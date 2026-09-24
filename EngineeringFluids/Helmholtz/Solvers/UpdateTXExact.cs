using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using System;

public static partial class Update
{
    // Exact (EOS-based) counterpart to UpdateTXFast.cs: instead of reading rhoL/rhoV/Psat off
    // the pre-fitted ancillary polynomials, this solves the real phase-equilibrium conditions
    // (equal pressure, equal fugacity) against the actual Helmholtz EOS via SolveAtTFast. Not
    // rejected within the last ~1 K of Tc or ~15-20 K of the triple point the way the ancillary
    // version's reliable range is bounded (see UpdateTXFast.cs remarks) - accuracy here comes
    // from Newton convergence against the EOS itself, not curve-fit quality, so it stays
    // reliable much closer to both endpoints (SolveAtTFast's own convergence checks are the
    // natural failure signal, not an arbitrary cutoff). The tradeoff is cost: a 2-equation
    // damped Newton solve per call instead of three O(1) polynomial reads - use UpdateTX for
    // the hot path and this when correctness right up to the boundary matters, or as a
    // reference to validate/regenerate ancillary fits (including for a future fluid that has
    // no fits yet at all).
    //
    // Within the last ~0.5 K of Tc, CoolProp's own saturation solver either disagrees with
    // this EOS by a small amount or fails outright (confirmed: SharpFluids returns a null
    // Pressure there) - not a bug on either side, just the well-known difficulty of resolving
    // near-critical saturation properties precisely, however good the solver. SolveAtTFast's
    // own pressure/fugacity residuals stay excellent (~1e-8 Pa / ~1e-14) right up to 0.001 K
    // from Tc regardless - there just isn't an independent oracle left to confirm it against
    // that close.
    public static void UpdateTXExact(this AmmoniaDouble local, double tTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        double Tc = local.Critical.Temperature;
        double Ttriple = local.TripleLiquid.Temperature;

        if (!double.IsFinite(tTarget) || tTarget >= Tc)
            throw new InvalidOperationException($"UpdateTXExact invalid at/above Tc. T={tTarget} K.");
        if (tTarget <= Ttriple)
            throw new InvalidOperationException($"UpdateTXExact invalid at/below triple temperature. T={tTarget} K.");

        var sat = local.SolveAtTFast(tTarget);
        local.SetTwoPhase(sat, quality);
    }

    // Alias to match SharpFluids naming (Quality, Temperature) - see UpdateTXFast.cs's UpdateXT.
    public static void UpdateXTExact(this AmmoniaDouble local, double quality, double tTarget)
        => UpdateTXExact(local, tTarget, quality);
}
