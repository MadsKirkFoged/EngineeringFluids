using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Temperature + entropy. Below Tc the exact saturation state at T decides first: an entropy between the
    // saturated liquid and vapor entropies is two-phase. Otherwise s falls monotonically with density along the
    // isotherm ((ds/drho)_T = -(dp/dT)_rho / rho^2), so the density is found with Brent's method in ln(rho) on the
    // right branch: below the saturated-vapor density for gas, above the saturated-liquid density for liquid, and
    // the whole range above Tc.
    //
    // Labels: CoolProp's plain single-phase rule (its TS flash labels are recalculated with it afterwards). An entropy
    // within DomeEdgeBand of a saturated entropy counts as saturated. Exactly at Tc CoolProp has no TS flash, and
    // neither has this.
    public static void UpdateTS(this Ammonia local, Temperature temperature, SpecificEntropy entropy)
    {
        if (!entropy.HasValue())
            throw new ArgumentOutOfRangeException(nameof(entropy), "UpdateTS: entropy must be finite.");
        if (!double.IsFinite(temperature.Kelvin) || temperature < local.TripleLiquid.Temperature)
            throw new ArgumentOutOfRangeException(nameof(temperature),
                $"UpdateTS: temperature below the triple point ({local.TripleLiquid.Temperature.Kelvin} K) is not supported. T={temperature.Kelvin} K.");
        if (temperature == local.Critical.Temperature)
            throw new InvalidOperationException("UpdateTS: not defined exactly at the critical temperature (as in CoolProp).");

        local.ClearTwoPhase();
        SetTS(local, temperature, entropy);
    }

    // The TS flash itself, also the inner evaluator of the HS flash
    private static void SetTS(Ammonia local, Temperature T, SpecificEntropy entropy)
    {
        MolarMass M = local.MolarMass;
        Molarity rhoc = local.Critical.MolarDensity;
        double sTarget = entropy.JoulePerKilogramKelvin;

        // Densities searched: from far below any gas to well above the densest liquid on the EOS's range
        double lnRhoMin = Math.Log(1e-10), lnRhoMax = Math.Log((6.0 * rhoc).MolesPerCubicMeter);

        if (T < local.Critical.Temperature && TrySaturationAtT(local, T) is SaturationSolver.SatResult sat)
        {
            double sL = new Ammonia { Temperature = T, Density = sat.RhomolarL * M }.Entropy.JoulePerKilogramKelvin;
            double sV = new Ammonia { Temperature = T, Density = sat.RhomolarV * M }.Entropy.JoulePerKilogramKelvin;
            double band = DomeEdgeBand * (sV - sL);

            if (sTarget >= sL - band && sTarget <= sV + band)
            {
                // Within the band of a saturated entropy the state is on that saturation line, q exactly 0 or 1. Inside
                // the HS flash the outer solve otherwise leaves noise-sized qualities (2.6e-12 seen on the 9 kPa liquid
                // line), and at low pressure the mixture density - and the compressibility CoolProp reports for it - is
                // sensitive to q even at that size.
                double q = sTarget <= sL + band ? 0.0 : sTarget >= sV - band ? 1.0 : (sTarget - sL) / (sV - sL);
                local.SetTwoPhase(sat, q);
                return;
            }

            if (sTarget > sV)
                lnRhoMax = Math.Log(sat.RhomolarV.MolesPerCubicMeter);
            else
                lnRhoMin = Math.Log(sat.RhomolarL.MolesPerCubicMeter);
        }

        // Single-phase from here on. Inside the HS flash the previous trial may have left this object two-phase,
        // and a two-phase object reports the mixture entropy regardless of the density being tried.
        local.ClearTwoPhase();

        double sScale = Math.Max(Math.Abs(sTarget), 1.0);
        double Residual(double lnRho)
        {
            local.Temperature = T;
            local.Density = Molarity.FromMolesPerCubicMeter(Math.Exp(lnRho)) * M;
            return (local.Entropy.JoulePerKilogramKelvin - sTarget) / sScale;
        }

        double flo = Residual(lnRhoMin), fhi = Residual(lnRhoMax);
        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"UpdateTS: could not bracket a density at T={T.Kelvin} K, s={sTarget} J/kg/K (outside the supported range).");

        double root = RootFinder.Brent(Residual, lnRhoMin, lnRhoMax, flo, fhi, xTol: 1e-14);
        if (!(Math.Abs(Residual(root)) <= 1e-9))
            throw new InvalidOperationException($"UpdateTS: did not converge. T={T.Kelvin} K, s={sTarget} J/kg/K.");
    }
}
