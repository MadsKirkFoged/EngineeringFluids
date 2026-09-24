using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;

public static partial class Update
{
    // Fast (AmmoniaDouble) pressure + vapor-quality flash: sets a two-phase saturated state
    // at pTarget with mass quality `quality` (0 = saturated liquid, 1 = saturated vapor).
    //
    // Unlike UpdatePT's Fast solver, this does NOT run an EOS-based Newton refinement: it
    // reads the saturation temperature and the two branch densities straight from the fast
    // ancillary polynomials (SaturationTemperatureFast / BubbleDensityFast / DewDensityFast).
    // Measured against CoolProp (via SharpFluids) across the full pressure range (from just
    // above the triple-point pressure up to within ~1 K of Tc), this already matches to
    // ~1e-6..1e-7 relative on Tsat/rhoL/rhoV - because the ancillaries were fit directly
    // against the same saturation curve this EOS is meant to reproduce, a coupled two-density
    // Newton solve (like the slow Ammonia.UpdatePX's SaturationSolver.SolveAtP) buys no
    // measurable accuracy here for a lot more per-call cost. Downstream Enthalpy/Entropy/
    // InternalEnergy are then computed by AmmoniaDouble's existing SetTwoPhase/SatLiquidState/
    // SatVaporState machinery, which evaluates the full EOS at (Tsat, rhoL) and (Tsat, rhoV) -
    // so their accuracy also inherits directly from how well the ancillaries track the true
    // saturation curve, not from any separate approximation. [benchmark-guided]
    //
    // The one place this genuinely breaks down is the last ~1 K below Tc, where rhoL and rhoV
    // converge onto the critical density and any small ancillary error is amplified in
    // relative terms - the exact same critical-point fragility already found and documented
    // for UpdatePT (both CoolProp's own saturation solver and this EOS struggle there). That
    // region is explicitly rejected below rather than silently returning a bad answer.
    public static void UpdatePX(this AmmoniaDouble local, double pTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new System.ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        double Tc = local.Critical.Temperature;
        double Pc = local.Critical.Pressure;
        double Ptriple = local.TripleLiquid.Pressure;

        if (!double.IsFinite(pTarget) || pTarget < Ptriple)
            throw new System.ArgumentOutOfRangeException(nameof(pTarget),
                $"UpdatePX: pressure below the triple-point pressure ({Ptriple} Pa) has no liquid-vapor equilibrium. P={pTarget} Pa.");

        if (pTarget >= Pc)
            throw new System.InvalidOperationException(
                $"UpdatePX invalid at/above the critical pressure. P={pTarget} Pa, Pc={Pc} Pa.");

        double Tsat = SaturationTemperatureFast.Temperature((float)pTarget);

        // Ancillary-based rhoL/rhoV lose accuracy fast in the last ~1 K below Tc (see remarks
        // above) - reject rather than silently return an inaccurate two-phase state there.
        if (Tc - Tsat < 1.0)
            throw new System.InvalidOperationException(
                $"UpdatePX: pressure too close to the critical point for a reliable fast ancillary-based " +
                $"two-phase solve (Tc-Tsat={Tc - Tsat:G3} K). P={pTarget} Pa. Use UpdatePT with a phase hint instead.");

        double rhoL = BubbleDensityFast.Density((float)Tsat);
        double rhoV = DewDensityFast.Density((float)Tsat);

        var sat = new SaturationSolver.SatResultDouble(Tsat, pTarget, rhoL, rhoV);
        local.SetTwoPhase(sat, quality);
    }
}
