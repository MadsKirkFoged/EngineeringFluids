using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringUnits;
using System;

namespace EngineeringFluids.Helmholtz.Solvers;

public static class SaturationSolver
{
    public sealed record SatResult(
        Temperature T,
        Pressure Psat,
        Molarity RhomolarL,
        Molarity RhomolarV);

    public static SatResult SolveAtT(this Ammonia local,Temperature T)
    {
        if (T == null)
            throw new ArgumentNullException(nameof(T));

        if (T <= local.TripleLiquid.Temperature || T >= local.Critical.Temperature)
            throw new ArgumentOutOfRangeException(nameof(T), "T must be between triple and critical for saturation.");

        var rhoL0 = LiquidDensity.CalculateDensity(T);
        var rhoV0 = VaporDensity.CalculateDensity(T);

        if (rhoL0 == null || rhoV0 == null)
            throw new InvalidOperationException("Ancillary density returned null; check temperature bounds.");

        double rhoL = rhoL0.MolesPerCubicMeter;
        double rhoV = rhoV0.MolesPerCubicMeter;

        if (rhoV <= 0 || rhoL <= 0 || rhoV >= rhoL)
        {
            rhoV = Math.Max(1e-8, rhoV);
            rhoL = Math.Max(rhoV * 10.0, rhoL);
        }

        const int maxIts = 80;
        const double tolP = 1e-2;
        const double tolLnPhi = 1e-12;

        Ammonia State(Temperature t, double rhomolar)
        {
            return new Ammonia
            {
                Temperature = t,
                Density = Molarity.FromMolesPerCubicMeter(rhomolar) * local.MolarMass
            };
        }

        // log variables: rhoV=exp(u), rhoL=exp(u+w)
        double u = Math.Log(rhoV);
        double w = Math.Log(rhoL / rhoV);

        double rhoRed = local.Critical.MolarDensity.MolesPerCubicMeter;
        double rhoMin = 1e-12;
        double rhoMax = 10.0 * rhoRed;

        for (int iter = 0; iter < maxIts; iter++)
        {
            rhoV = Math.Exp(u);
            rhoL = Math.Exp(u + w);

            var V = State(T, rhoV);
            var L = State(T, rhoL);

            double pV = V.Pressure.Pascal;
            double pL = L.Pressure.Pascal;

            double lnphiV = V.LNFugacityCoefficient;
            double lnphiL = L.LNFugacityCoefficient;

            double F1 = pL - pV;
            double F2 = lnphiL - lnphiV;

            if (Math.Abs(F1) < tolP && Math.Abs(F2) < tolLnPhi)
            {
                var psat = Pressure.FromPascal(0.5 * (pL + pV));
                return new SatResult(
                    T,
                    psat,
                    Molarity.FromMolesPerCubicMeter(rhoL),
                    Molarity.FromMolesPerCubicMeter(rhoV)
                );
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


    public static SatResult SolveAtP(this Ammonia local, Pressure P)
    {
        if (P == null)
            throw new ArgumentNullException(nameof(P));

        Temperature Tguess = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(P);

        double Ttr = local.TripleLiquid.Temperature.Kelvin;
        double Tc = local.Critical.Temperature.Kelvin;

        double T0 = Math.Clamp(Tguess.Kelvin, Ttr + 1e-3, Tc - 1e-3);

        double span = 5.0;
        double Tlo = Math.Max(Ttr + 1e-3, T0 - span);
        double Thi = Math.Min(Tc - 1e-3, T0 + span);

        SatResult satLo = local.SolveAtT(Temperature.FromKelvin(Tlo));
        SatResult satHi = local.SolveAtT(Temperature.FromKelvin(Thi));

        double fLo = satLo.Psat.Pascal - P.Pascal;
        double fHi = satHi.Psat.Pascal - P.Pascal;

        for (int expand = 0; expand < 40 && Math.Sign(fLo) == Math.Sign(fHi); expand++)
        {
            span *= 1.6;

            Tlo = Math.Max(Ttr + 1e-3, T0 - span);
            Thi = Math.Min(Tc - 1e-3, T0 + span);

            satLo = local.SolveAtT(Temperature.FromKelvin(Tlo));
            satHi = local.SolveAtT(Temperature.FromKelvin(Thi));

            fLo = satLo.Psat.Pascal - P.Pascal;
            fHi = satHi.Psat.Pascal - P.Pascal;

            if (Tlo <= Ttr + 1e-3 && Thi >= Tc - 1e-3 && Math.Sign(fLo) == Math.Sign(fHi))
                throw new InvalidOperationException($"Could not bracket Tsat for P={P.Pascal} Pa.");
        }

        if (Math.Sign(fLo) == Math.Sign(fHi))
            throw new InvalidOperationException("Could not bracket saturation temperature for given pressure.");

        const int maxIts = 80;
        const double tolP = 1e-1;

        SatResult midSat = satLo;

        for (int i = 0; i < maxIts; i++)
        {
            double Tmid = 0.5 * (Tlo + Thi);
            midSat = local.SolveAtT(Temperature.FromKelvin(Tmid));

            double fMid = midSat.Psat.Pascal - P.Pascal;

            if (Math.Abs(fMid) < tolP)
                return midSat;

            if (Math.Sign(fMid) == Math.Sign(fLo))
            {
                Tlo = Tmid;
                fLo = fMid;
            }
            else
            {
                Thi = Tmid;
                fHi = fMid;
            }
        }

        return midSat;
    }

}