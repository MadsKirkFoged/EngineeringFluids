using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;

public static partial class Update
{
    // Exact (EOS-based) counterpart to UpdateTX.cs: instead of reading rhoL/rhoV/Psat off
    // the ancillary tables, this solves the real phase-equilibrium conditions (equal pressure,
    // equal fugacity) against the actual Helmholtz EOS via SolveAtT, and is not rejected in the
    // last millikelvin below Tc the way UpdateTX is - SolveAtT converges to within 6 microK of Tc
    // (its own checks are the failure signal, not a cutoff). The tradeoff is cost: a 2-equation
    // damped Newton solve per call instead of three O(1) table reads - use UpdateTX for the hot
    // path and this for SharpFluids compatibility, or as a reference to validate/regenerate the
    // ancillary tables (including for a future fluid that has none yet).
    //
    // Above ~405.13 K CoolProp's own saturation solver goes wrong (non-monotonic, then garbage
    // at 405.17 K and a native crash), so there is no CoolProp reference that close; the Gao
    // dome itself continues to Tc (checked by continuation, see the ancillary tables' remarks).
    public static void UpdateTXExact(this Ammonia local, Temperature tTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        Temperature Tc = local.Critical.Temperature;
        Temperature Ttriple = local.TripleLiquid.Temperature;

        if (!double.IsFinite(tTarget.Kelvin) || tTarget >= Tc)
            throw new InvalidOperationException($"UpdateTXExact invalid at/above Tc. T={tTarget.Kelvin} K.");
        if (tTarget < Ttriple)
            throw new InvalidOperationException($"UpdateTXExact invalid below triple temperature. T={tTarget.Kelvin} K.");

        var sat = local.SolveAtT(tTarget);
        local.SetTwoPhase(sat, quality);
    }

    // Alias to match SharpFluids naming (Quality, Temperature) - see UpdateTX.cs's UpdateXT.
    public static void UpdateXTExact(this Ammonia local, double quality, Temperature tTarget)
        => UpdateTXExact(local, tTarget, quality);
}
