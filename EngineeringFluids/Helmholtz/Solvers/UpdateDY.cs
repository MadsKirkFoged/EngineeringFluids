using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;

public static partial class Update
{
    // Density + enthalpy and density + entropy. Along an isochore both h and s rise monotonically with T - through
    // the dome too (there the DT flash gives the lever-rule mixture) - so T is found with Brent's method over
    // DT flashes, between the triple point and CoolProp's search limit. Labels: CoolProp's plain single-phase rule
    // (its HSU_D flash ends with recalculate_singlephase_phase), Twophase inside the dome.
    public static void UpdateDH(this Ammonia local, Density density, Enthalpy enthalpy)
    {
        if (!enthalpy.HasValue())
            throw new ArgumentOutOfRangeException(nameof(enthalpy), "UpdateDH: enthalpy must be finite.");
        UpdateDY(local, density, enthalpy.JoulePerKilogram, FlashProperty.Enthalpy, nameof(UpdateDH));
    }

    public static void UpdateDS(this Ammonia local, Density density, SpecificEntropy entropy)
    {
        if (!entropy.HasValue())
            throw new ArgumentOutOfRangeException(nameof(entropy), "UpdateDS: entropy must be finite.");
        UpdateDY(local, density, entropy.JoulePerKilogramKelvin, FlashProperty.Entropy, nameof(UpdateDS));
    }

    private static void UpdateDY(Ammonia local, Density density, double yTarget, FlashProperty y, string name)
    {
        if (!double.IsFinite(density.KilogramPerCubicMeter) || density <= Density.Zero)
            throw new ArgumentOutOfRangeException(nameof(density), $"{name}: density must be positive. D={density.KilogramPerCubicMeter} kg/m3.");

        local.ClearTwoPhase();
        double yScale = Math.Max(Math.Abs(yTarget), 1.0);

        double Residual(double TK)
        {
            SetDT(local, density, Temperature.FromKelvin(TK));
            return (Read(local, y) - yTarget) / yScale;
        }

        double lo = local.TripleLiquid.Temperature.Kelvin, hi = ExactFlashTmax.Kelvin;
        double flo = Residual(lo), fhi = Residual(hi);
        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"{name}: could not bracket a solution in [{lo:G6} K, {hi:G6} K]. D={density.KilogramPerCubicMeter} kg/m3, {Describe(y, yTarget)} is likely outside the supported range.");

        double root = RootFinder.Brent(Residual, lo, hi, flo, fhi, xTol: 1e-10);
        if (!(Math.Abs(Residual(root)) <= 1e-9))
            throw new InvalidOperationException($"{name}: did not converge. D={density.KilogramPerCubicMeter} kg/m3, {Describe(y, yTarget)}.");

        if (local.Quality < 0)
            local.SetPhase(local.SinglePhaseLabel());
    }
}
