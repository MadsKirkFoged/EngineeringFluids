using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;

public static partial class Update
{
    // Temperature + vapor-quality flash: sets a two-phase saturated state
    // at tTarget with mass quality `quality` (0 = saturated liquid, 1 = saturated vapor).
    //
    // Even cheaper than UpdatePX: T is already the independent variable of the ancillary
    // tables, so this goes straight to Psat(T)/rhoL(T)/rhoV(T) - three O(1) table reads, no
    // iteration of any kind. Measured against UpdateTXExact (p, rho, h, s, u): <= 5e-10 relative
    // from the triple point to 10 mK below Tc, 2e-8 between 10 mK and 1 mK; the last millikelvin
    // is rejected below, as in UpdatePX.
    public static void UpdateTX(this Ammonia local, Temperature tTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        Temperature Tc = local.Critical.Temperature;
        Temperature Ttriple = local.TripleLiquid.Temperature;

        if (!double.IsFinite(tTarget.Kelvin) || tTarget >= Tc)
            throw new InvalidOperationException($"UpdateTX invalid at/above Tc. T={tTarget.Kelvin} K.");

        if (tTarget < Ttriple)
            throw new InvalidOperationException($"UpdateTX invalid below triple temperature. T={tTarget.Kelvin} K.");

        // Within the last millikelvin below Tc even exact rhoL/rhoV are ill-conditioned (see remarks
        // above) - reject rather than silently return an inaccurate two-phase state there.
        if (Tc - tTarget < CriticalMargin)
            throw new InvalidOperationException(
                $"UpdateTX: temperature too close to the critical point for a reliable fast ancillary-based " +
                $"two-phase solve (Tc-T={(Tc - tTarget).Kelvin:G3} K). T={tTarget.Kelvin} K. Use UpdatePT with a phase hint instead.");

        Pressure psat = SaturationPressure.Pressure(tTarget);
        Molarity rhoL = BubbleDensity.Density(tTarget);
        Molarity rhoV = DewDensity.Density(tTarget);

        var sat = new SaturationSolver.SatResult(tTarget, psat, rhoL, rhoV);
        local.SetTwoPhase(sat, quality);
    }

    // Alias to match SharpFluids naming (Quality, Temperature).
    public static void UpdateXT(this Ammonia local, double quality, Temperature tTarget)
        => UpdateTX(local, tTarget, quality);
}
