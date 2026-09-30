using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;

public static partial class Update
{
    // Enthalpy + entropy, the way CoolProp's HS flash does it: along an isentrope h rises monotonically with T
    // (dh = v dp there, and p rises with T), so T is found with Brent's method over TS flashes between the triple
    // point and the EOS's maximum temperature. Like CoolProp, inputs that would need a higher temperature fail.
    // Labels come from the final TS flash.
    public static void UpdateHS(this Ammonia local, Enthalpy enthalpy, SpecificEntropy entropy)
    {
        if (!enthalpy.HasValue())
            throw new ArgumentOutOfRangeException(nameof(enthalpy), "UpdateHS: enthalpy must be finite.");
        if (!entropy.HasValue())
            throw new ArgumentOutOfRangeException(nameof(entropy), "UpdateHS: entropy must be finite.");

        local.ClearTwoPhase();

        double hTarget = enthalpy.JoulePerKilogram;
        double hScale = Math.Max(Math.Abs(hTarget), 1.0);
        double Tc = local.Critical.Temperature.Kelvin;

        double Residual(double TK)
        {
            // The TS flash is undefined exactly at Tc; the solution is continuous there, so step off it
            if (TK == Tc)
                TK = Math.BitIncrement(TK);
            SetTS(local, Temperature.FromKelvin(TK), entropy);
            return (local.Enthalpy.JoulePerKilogram - hTarget) / hScale;
        }

        bool TryResidual(double TK, out double f)
        {
            try { f = Residual(TK); return true; }
            catch (InvalidOperationException) { f = double.NaN; return false; }
        }

        double lo = local.TripleLiquid.Temperature.Kelvin, hi = local.LimitTemperatureMax.Kelvin;
        if (!TryResidual(lo, out double flo))
            throw new InvalidOperationException(
                $"UpdateHS: no state with s={entropy.JoulePerKilogramKelvin} J/kg/K at the triple-point temperature; the input is outside the supported range.");

        // A low entropy at high T would need a density beyond the EOS's range, so the TS flash fails there. Where it
        // works is a contiguous range from the triple point up (lower T needs lower density for the same s), so the
        // upper end of the bracket is pulled down to the highest temperature that still has a TS state.
        if (!TryResidual(hi, out double fhi))
        {
            double good = lo, bad = hi;
            for (int i = 0; i < 60; i++)
            {
                double mid = 0.5 * (good + bad);
                if (TryResidual(mid, out _))
                    good = mid;
                else
                    bad = mid;
            }
            hi = good;
            fhi = Residual(hi);
        }
        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"UpdateHS: could not bracket a solution in [{lo:G6} K, {hi:G6} K]. h={hTarget} J/kg, s={entropy.JoulePerKilogramKelvin} J/kg/K is outside the supported range.");

        double root = RootFinder.Brent(Residual, lo, hi, flo, fhi, xTol: 1e-10);
        if (!(Math.Abs(Residual(root)) <= 1e-9))
            throw new InvalidOperationException($"UpdateHS: did not converge. h={hTarget} J/kg, s={entropy.JoulePerKilogramKelvin} J/kg/K.");
    }
}
