using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using System;

public static partial class Update
{
    // Fast (AmmoniaDouble) temperature + vapor-quality flash: sets a two-phase saturated state
    // at tTarget with mass quality `quality` (0 = saturated liquid, 1 = saturated vapor).
    //
    // Even cheaper than UpdatePXFast: T is already the independent variable the ancillaries
    // are fit against, so this skips UpdatePX's own P->T bisection (SaturationTemperatureFast)
    // entirely and goes straight to Psat(T)/rhoL(T)/rhoV(T) - three direct O(1) ancillary
    // polynomial reads, no iteration of any kind. Same accuracy profile and near-critical
    // fragility as UpdatePXFast (rhoL/rhoV both converge onto the critical density in the last
    // ~1 K below Tc, amplifying any small ancillary error in relative terms) - rejected below
    // for the same reason. [benchmark-guided]
    //
    // The reduced-accuracy band is wider than that hard 1 K cutoff, though: sampling T
    // directly (rather than solving it from P, as UpdatePX's own sweep does) showed
    // rhoL/rhoV relative error up to ~3e-4 as far as ~10-11 K below Tc, and Psat's own fit
    // breaching 1e-5 as far as ~15 K above the triple point - the same ancillary-fit dip
    // already known from UpdatePX's pressure-domain sweep (documented there as the
    // "0.8-0.95 Pc dip"), just directly exposed here across a wider temperature margin than
    // that P-domain sampling ever exercised. Not rejected outright (unlike the 1 K-from-Tc
    // band, where the ancillaries are unconditionally unreliable) since accuracy degrades
    // gradually rather than breaking down - callers needing tight tolerances within ~15 K of
    // either endpoint should verify against UpdatePT instead. [benchmark-guided]
    public static void UpdateTX(this AmmoniaDouble local, double tTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        double Tc = local.Critical.Temperature;
        double Ttriple = local.TripleLiquid.Temperature;

        if (!double.IsFinite(tTarget) || tTarget >= Tc)
            throw new InvalidOperationException($"UpdateTX invalid at/above Tc. T={tTarget} K.");

        if (tTarget < Ttriple)
            throw new InvalidOperationException($"UpdateTX invalid below triple temperature. T={tTarget} K.");

        // Ancillary-based rhoL/rhoV lose accuracy fast in the last ~1 K below Tc (see remarks
        // above) - reject rather than silently return an inaccurate two-phase state there.
        if (Tc - tTarget < 1.0)
            throw new InvalidOperationException(
                $"UpdateTX: temperature too close to the critical point for a reliable fast ancillary-based " +
                $"two-phase solve (Tc-T={Tc - tTarget:G3} K). T={tTarget} K. Use UpdatePT with a phase hint instead.");

        double psat = SaturationPressureFast.Pressure((float)tTarget);
        double rhoL = BubbleDensityFast.Density((float)tTarget);
        double rhoV = DewDensityFast.Density((float)tTarget);

        var sat = new SaturationSolver.SatResultDouble(tTarget, psat, rhoL, rhoV);
        local.SetTwoPhase(sat, quality);
    }

    // Alias to match SharpFluids naming (Quality, Temperature).
    public static void UpdateXT(this AmmoniaDouble local, double quality, double tTarget)
        => UpdateTX(local, tTarget, quality);
}
