using EngineeringFluids.Fluids;
using System;

namespace EngineeringFluids.Helmholtz.Solvers;

public static class SaturationSolver
{
    public sealed record SatResult(
    double T,
    double Psat,
    double RhomolarL,
    double RhomolarV);

    // Rigorous saturation solve: equal pressure + equal fugacity between the two branches,
    // damped Newton in log-density space u=ln(rhoV), w=ln(rhoL/rhoV). Unlike
    // UpdateTX/UpdatePX, this does NOT read rhoL/rhoV/Psat off a
    // pre-fitted ancillary curve - it solves the actual phase-equilibrium conditions against
    // the real EOS every call, so its accuracy is bounded only by EOS/Newton convergence
    // (~1e-10 relative), not by how well a polynomial happens to track the true saturation
    // line. That's what makes it valuable beyond just "a slower UpdateTX": it's the same kind
    // of solver that would be needed to GENERATE ancillary fits for any new fluid in the first
    // place, since there is nothing else to fit those curves against.
    //
    // The one piece that stays fluid-specific is the initial guess: BubbleDensity/
    // DewDensity give a good starting point for ammonia because those fits already exist.
    // A brand-new fluid without ancillary fits yet would need a different bootstrap (e.g. a
    // corresponding-states or Clausius-Clapeyron estimate, or continuation from a temperature
    // where a guess IS available) - only this seed step, not the Newton iteration itself,
    // would need to change per fluid.
    public static SatResult SolveAtT(this Ammonia local, double T)
    {
        double Ttriple = local.TripleLiquid.Temperature;
        double Tc = local.Critical.Temperature;

        if (!double.IsFinite(T) || T <= Ttriple || T >= Tc)
            throw new ArgumentOutOfRangeException(nameof(T), "T must be between triple and critical for saturation.");

        double M = local.MolarMass;
        double rhoRed = local.Critical.MolarDensity;

        double rhoL = BubbleDensity.Density((float)T);
        double rhoV = DewDensity.Density((float)T);

        if (rhoV <= 0 || rhoL <= 0 || rhoV >= rhoL)
        {
            rhoV = Math.Max(1e-8, rhoV);
            rhoL = Math.Max(rhoV * 10.0, rhoL);
        }

        const int maxIts = 80;
        const double tolP = 1e-2;
        // 1e-12 turned out unreachable at some low-T points: the
        // residual settles at the floating-point noise floor (observed oscillating around
        // ~1e-12 to 1e-6, never both below 1e-12 at once) once Newton has genuinely converged
        // - du/dw were already ~1e-15, i.e. pure rounding noise, not real correction. 1e-10 is
        // still far tighter than any practical need (ln(phi) resolved to 10 digits) and clears
        // that floor everywhere tested. [benchmark-guided]
        const double tolLnPhi = 1e-10;

        static Ammonia State(double t, double rhomolar, double m) => new Ammonia
        {
            Temperature = t,
            Density = rhomolar * m
        };

        double u = Math.Log(rhoV);
        double w = Math.Log(rhoL / rhoV);

        double rhoMin = 1e-12;
        double rhoMax = 10.0 * rhoRed;

        for (int iter = 0; iter < maxIts; iter++)
        {
            rhoV = Math.Exp(u);
            rhoL = Math.Exp(u + w);

            var V = State(T, rhoV, M);
            var L = State(T, rhoL, M);

            double pV = V.Pressure;
            double pL = L.Pressure;

            double lnphiV = V.LNFugacityCoefficient;
            double lnphiL = L.LNFugacityCoefficient;

            double F1 = pL - pV;
            double F2 = lnphiL - lnphiV;

            if (Math.Abs(F1) < tolP && Math.Abs(F2) < tolLnPhi)
            {
                double psat = 0.5 * (pL + pV);
                return new SatResult(T, psat, rhoL, rhoV);
            }

            double dpL = L.dp_drhomolar_constT_SI;
            double dpV = V.dp_drhomolar_constT_SI;

            double dlnphiL = L.dLnPhi_dRhomolar_constT_SI;
            double dlnphiV = V.dLnPhi_dRhomolar_constT_SI;

            double a11 = dpL * rhoL - dpV * rhoV;
            double a12 = dpL * rhoL;
            double a21 = dlnphiL * rhoL - dlnphiV * rhoV;
            double a22 = dlnphiL * rhoL;

            double det = a11 * a22 - a12 * a21;
            if (!double.IsFinite(det) || Math.Abs(det) < 1e-30)
                throw new InvalidOperationException($"Saturation Jacobian singular at iter {iter}.");

            double b1 = -F1;
            double b2 = -F2;

            double du = (b1 * a22 - a12 * b2) / det;
            double dw = (a11 * b2 - b1 * a21) / det;

            // The noise floor above varies from point to point (observed ~1e-12 at one T,
            // ~1.7e-10 at another) - no fixed tolLnPhi clears it everywhere. Once Newton's own
            // proposed correction is already negligible (~1e-15, pure rounding noise) while the
            // pressure residual is genuinely converged, further iteration can only bounce
            // between floating-point-adjacent values, never actually improving F2 - accept the
            // current point rather than exhaust maxIts chasing a residual below the noise
            // floor. [benchmark-guided]
            if (Math.Abs(F1) < tolP && Math.Abs(du) < 1e-11 && Math.Abs(dw) < 1e-11)
            {
                double psat = 0.5 * (pL + pV);
                return new SatResult(T, psat, rhoL, rhoV);
            }

            double lambda = 1.0;
            bool accepted = false;

            for (int damp = 0; damp < 30; damp++)
            {
                double uTry = u + lambda * du;
                double wTry = w + lambda * dw;

                double rhoVTry = Math.Exp(uTry);
                double rhoLTry = Math.Exp(uTry + wTry);

                bool ok =
                    rhoVTry > rhoMin &&
                    rhoLTry > rhoVTry &&
                    rhoVTry < rhoMax &&
                    rhoLTry < rhoMax &&
                    wTry > 0.0;

                if (ok)
                {
                    u = uTry;
                    w = wTry;
                    accepted = true;
                    break;
                }

                lambda *= 0.5;
            }

            if (!accepted)
                throw new InvalidOperationException($"Saturation damping failed at iter {iter}. u={u}, w={w}");
        }

        throw new InvalidOperationException("Saturation solver did not converge.");
    }
}