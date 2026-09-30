using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using System;

namespace EngineeringFluids.Helmholtz.Solvers;

public static class SaturationSolver
{
    public sealed record SatResult(
    Temperature T,
    Pressure Psat,
    Molarity RhomolarL,
    Molarity RhomolarV);

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
    public static SatResult SolveAtT(this Ammonia local, Temperature T)
    {
        Temperature Ttriple = local.TripleLiquid.Temperature;
        Temperature Tc = local.Critical.Temperature;

        if (!double.IsFinite(T.Kelvin) || T < Ttriple || T >= Tc)
            throw new ArgumentOutOfRangeException(nameof(T), "T must be between triple and critical for saturation.");

        MolarMass M = local.MolarMass;
        Molarity rhoRed = local.Critical.MolarDensity;

        Molarity rhoL = BubbleDensity.Density(T);
        Molarity rhoV = DewDensity.Density(T);

        // ln(rhoL/rhoV) of the ancillaries, the yardstick NonTrivial checks the converged dome width against (0 = no usable
        // ancillary width: both fits end at rhoc, so right at Tc they may meet or cross)
        double wAncillary = rhoV > Molarity.Zero && rhoL > rhoV ? Math.Log((double)(rhoL / rhoV)) : 0.0;

        if (rhoV <= Molarity.Zero || rhoL <= Molarity.Zero || rhoV >= rhoL)
        {
            rhoV = Molarity.Max(Molarity.FromMolesPerCubicMeter(1e-8), rhoV);
            rhoL = Molarity.Max(rhoV * 10.0, rhoL);
        }

        const int maxIts = 80;
        // Relative, not absolute: the former fixed 1e-2 Pa is 1.6e-6 of the triple-point pressure, which capped
        // SolveAtP's accuracy there. Near the triple point the residual still jitters between ~1e-11 and ~1e-8 from
        // one iterate to the next (measured), so 1e-10 was not reliably reachable; 1e-9 is, and equals the old
        // tolerance at ~10 MPa.
        const double tolPRel = 1e-9;
        // 1e-12 turned out unreachable at some low-T points: the
        // residual settles at the floating-point noise floor (observed oscillating around
        // ~1e-12 to 1e-6, never both below 1e-12 at once) once Newton has genuinely converged
        // - du/dw were already ~1e-15, i.e. pure rounding noise, not real correction. 1e-10 is
        // still far tighter than any practical need (ln(phi) resolved to 10 digits) and clears
        // that floor everywhere tested. [benchmark-guided]
        const double tolLnPhi = 1e-10;
        // ln(rhoL/rhoV) below which the dome counts as narrow (rhoL/rhoV < 1.105, T within ~9 mK of Tc)
        const double NarrowDome = 0.1;

        static Ammonia State(Temperature t, Molarity rhomolar, MolarMass m) => new Ammonia
        {
            Temperature = t,
            Density = rhomolar * m
        };

        // Newton runs in log-density space, so the logarithms are taken of the plain mol/m3 values.
        double u = Math.Log(rhoV.MolesPerCubicMeter);
        double w = Math.Log((double)(rhoL / rhoV));

        int narrowPolish = 0;
        Molarity rhoMin = Molarity.FromMolesPerCubicMeter(1e-12);
        Molarity rhoMax = 10.0 * rhoRed;

        for (int iter = 0; iter < maxIts; iter++)
        {
            rhoV = Molarity.FromMolesPerCubicMeter(Math.Exp(u));
            rhoL = Molarity.FromMolesPerCubicMeter(Math.Exp(u + w));

            var V = State(T, rhoV, M);
            var L = State(T, rhoL, M);

            Pressure pV = V.Pressure;
            Pressure pL = L.Pressure;

            double lnphiV = V.LNFugacityCoefficient;
            double lnphiL = L.LNFugacityCoefficient;

            Pressure F1 = pL - pV;
            double F2 = lnphiL - lnphiV;

            bool residualsConverged = F1.Abs() < tolPRel * pV && Math.Abs(F2) < tolLnPhi;
            if (residualsConverged && w > NarrowDome)
            {
                Pressure psat = 0.5 * (pL + pV);
                return NonTrivial(new SatResult(T, psat, rhoL, rhoV), wAncillary);
            }

            MolarEnergy dpL = L.dp_drhomolar_constT;
            MolarEnergy dpV = V.dp_drhomolar_constT;

            var dlnphiL = L.dLnPhi_dRhomolar_constT;
            var dlnphiV = V.dLnPhi_dRhomolar_constT;

            Pressure a11 = dpL * rhoL - dpV * rhoV;
            Pressure a12 = dpL * rhoL;
            double a21 = (double)(dlnphiL * rhoL - dlnphiV * rhoV);
            double a22 = (double)(dlnphiL * rhoL);

            Pressure det = a11 * a22 - a12 * a21;
            if (!det.HasValue() || det.Abs() < Pressure.FromPascal(1e-30))
                throw new InvalidOperationException($"Saturation Jacobian singular at iter {iter}.");

            Pressure b1 = -F1;
            double b2 = -F2;

            double du = (double)((b1 * a22 - a12 * b2) / det);
            double dw = (double)((a11 * b2 - b1 * a21) / det);

            // In a narrow dome (the last ~9 mK below Tc) p and ln(phi) are nearly flat across it, so converged residuals
            // no longer pin the densities: they still left rhoL, rhoV ~1e-5 off 0.4 mK below Tc. There Newton keeps
            // going until its width step is negligible - or, since that step bottoms out at the rounding noise of the
            // ill-conditioned problem closer to Tc, for at most three more (quadratically converging) steps.
            if (residualsConverged && (Math.Abs(dw) <= 1e-9 * w || ++narrowPolish > 3))
            {
                Pressure psat = 0.5 * (pL + pV);
                return NonTrivial(new SatResult(T, psat, rhoL, rhoV), wAncillary);
            }

            // The noise floor above varies from point to point (observed ~1e-12 at one T,
            // ~1.7e-10 at another) - no fixed tolLnPhi clears it everywhere. Once Newton's own
            // proposed correction is already negligible (~1e-15, pure rounding noise) while the
            // pressure residual is genuinely converged, further iteration can only bounce
            // between floating-point-adjacent values, never actually improving F2 - accept the
            // current point rather than exhaust maxIts chasing a residual below the noise
            // floor. [benchmark-guided]
            if (F1.Abs() < tolPRel * pV && Math.Abs(du) < 1e-11 && Math.Abs(dw) < 1e-11)
            {
                Pressure psat = 0.5 * (pL + pV);
                return NonTrivial(new SatResult(T, psat, rhoL, rhoV), wAncillary);
            }

            double lambda = 1.0;
            bool accepted = false;

            for (int damp = 0; damp < 30; damp++)
            {
                double uTry = u + lambda * du;
                double wTry = w + lambda * dw;

                Molarity rhoVTry = Molarity.FromMolesPerCubicMeter(Math.Exp(uTry));
                Molarity rhoLTry = Molarity.FromMolesPerCubicMeter(Math.Exp(uTry + wTry));

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

    // Equal pressure and equal fugacity are also satisfied trivially by rhoL == rhoV, and just below Tc the damped Newton
    // could slide onto it or onto a nearby spurious root. With the old float ancillaries as the start that happened from
    // 405.24 K to 405.38 K: rhoL/rhoV = 1.0002 where the real dome is 1.378 wide, p 1.2e-4 off - and a bare
    // rhoL > 1.000001 rhoV check let it through. The ancillary dome width ln(rhoL/rhoV) is now accurate to ~1e-4 of
    // itself even 10 microK below Tc, so a real solution lies within a factor 2 of it; anything else is a failure, not a
    // saturation state.
    private static SatResult NonTrivial(SatResult sat, double wAncillary)
    {
        double w = Math.Log((double)(sat.RhomolarL / sat.RhomolarV));
        if (!(sat.RhomolarL > 1.000001 * sat.RhomolarV) || (wAncillary > 0 && !(w > 0.5 * wAncillary && w < 2 * wAncillary)))
            throw new InvalidOperationException($"Saturation solver converged to a trivial or spurious root at T={sat.T.Kelvin} K: ln(rhoL/rhoV) = {w}, ancillary {wAncillary}.");
        return sat;
    }

    // Saturation state at a given pressure, solved against the EOS like SolveAtT (no ancillary accuracy limit):
    // Newton on T with the Clausius-Clapeyron slope dp/dT = (s_V - s_L) / (v_V - v_L) around SolveAtT, started
    // from the SaturationTemperature ancillary. That start is already within ~1e-8 K (2e-9 relative in p), so it
    // typically takes one or two SolveAtT calls. The returned Psat is the requested pressure itself.
    public static SatResult SolveAtP(this Ammonia local, Pressure p)
    {
        Pressure Ptriple = local.TripleLiquid.Pressure;
        Pressure Pc = local.Critical.Pressure;

        if (!double.IsFinite(p.Pascal) || p < Ptriple || p >= Pc)
            throw new ArgumentOutOfRangeException(nameof(p), "p must be between the triple-point and critical pressures for saturation.");

        Temperature Ttriple = local.TripleLiquid.Temperature;
        Temperature Tc = local.Critical.Temperature;
        Temperature lowest = Ttriple;
        Temperature highest = Tc - Temperature.FromKelvin(1e-9);
        MolarMass M = local.MolarMass;

        Temperature T = SaturationTemperature.Temperature(p).Clamp(lowest, highest);

        const int maxIts = 30;
        for (int iter = 0; iter < maxIts; iter++)
        {
            SatResult sat = local.SolveAtT(T);
            Pressure f = sat.Psat - p;

            // SolveAtT resolves Psat to ~1e-9 relative (its pressure tolerance), so asking for more only makes Newton
            // chase noise-sized steps (it did at ~8 kPa); 2e-9 still pins Tsat to ~1e-10 relative
            if (f.Abs() <= 2e-9 * p)
                return new SatResult(sat.T, p, sat.RhomolarL, sat.RhomolarV);

            var L = new Ammonia { Temperature = T, Density = sat.RhomolarL * M };
            var V = new Ammonia { Temperature = T, Density = sat.RhomolarV * M };
            MolarEntropy dS = V.MolarEntropy - L.MolarEntropy;
            UnknownUnit dV = 1.0 / sat.RhomolarV - 1.0 / sat.RhomolarL;

            Temperature step = f * dV / dS;
            T = (T - step).Clamp(lowest, highest);

            if (step.Abs() <= Temperature.FromKelvin(1e-10))
            {
                SatResult last = local.SolveAtT(T);
                if ((last.Psat - p).Abs() <= 1e-8 * p)
                    return new SatResult(last.T, p, last.RhomolarL, last.RhomolarV);
                break;
            }
        }

        throw new InvalidOperationException($"SolveAtP did not converge. P={p.Pascal} Pa.");
    }

    // Saturation state at T, the triple point included. Note that CoolProp's tabulated triple-point densities
    // (TripleLiquid/TripleVapor) do not belong to this EOS - at them it gives p = -2.6 MPa for the liquid - so the
    // triple point is solved like any other temperature (its EOS saturation pressure is 6055.8 Pa, vs the tabulated
    // 6091.2 Pa; CoolProp's own QT flash agrees with the EOS value).
    internal static SatResult SaturationAtT(this Ammonia local, Temperature T) => local.SolveAtT(T);

    // The critical point as a degenerate saturation state (both phases at the critical density).
    internal static SatResult CriticalPoint(Ammonia local)
        => new SatResult(local.Critical.Temperature, local.Critical.Pressure, local.Critical.MolarDensity, local.Critical.MolarDensity);
}