using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;

public static partial class Update
{
    // The ancillary-based two-phase flashes (UpdatePX/UpdateTX/UpdatePH/UpdatePS) are not
    // trusted within this distance below Tc - see the remarks on UpdatePX.
    private static readonly Temperature CriticalMargin = Temperature.FromKelvin(1.0);

    // Pressure + vapor-quality flash: sets a two-phase saturated state
    // at pTarget with mass quality `quality` (0 = saturated liquid, 1 = saturated vapor).
    //
    // Unlike UpdatePT, this does NOT run an EOS-based Newton refinement: it
    // reads the saturation temperature and the two branch densities straight from the fast
    // ancillary polynomials (SaturationTemperature / BubbleDensity / DewDensity).
    // Measured against CoolProp (via SharpFluids) across the full pressure range (from just
    // above the triple-point pressure up to within ~1 K of Tc), this already matches to
    // ~1e-6..1e-7 relative on Tsat/rhoL/rhoV - because the ancillaries were fit directly
    // against the same saturation curve this EOS is meant to reproduce, a coupled two-density
    // Newton solve (like SaturationSolver.SolveAtT) buys no
    // measurable accuracy here for a lot more per-call cost. Downstream Enthalpy/Entropy/
    // InternalEnergy are then computed by Ammonia's existing SetTwoPhase/SatLiquidState/
    // SatVaporState machinery, which evaluates the full EOS at (Tsat, rhoL) and (Tsat, rhoV) -
    // so their accuracy also inherits directly from how well the ancillaries track the true
    // saturation curve, not from any separate approximation. [benchmark-guided]
    //
    // The one place this genuinely breaks down is the last ~1 K below Tc, where rhoL and rhoV
    // converge onto the critical density and any small ancillary error is amplified in
    // relative terms - the exact same critical-point fragility already found and documented
    // for UpdatePT (both CoolProp's own saturation solver and this EOS struggle there). That
    // region is explicitly rejected below rather than silently returning a bad answer.
    public static void UpdatePX(this Ammonia local, Pressure pTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new System.ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        Temperature Tc = local.Critical.Temperature;
        Pressure Pc = local.Critical.Pressure;
        Pressure Ptriple = local.TripleLiquid.Pressure;

        if (!double.IsFinite(pTarget.Pascal) || pTarget < Ptriple)
            throw new System.ArgumentOutOfRangeException(nameof(pTarget),
                $"UpdatePX: pressure below the triple-point pressure ({Ptriple.Pascal} Pa) has no liquid-vapor equilibrium. P={pTarget.Pascal} Pa.");

        if (pTarget >= Pc)
            throw new System.InvalidOperationException(
                $"UpdatePX invalid at/above the critical pressure. P={pTarget.Pascal} Pa, Pc={Pc.Pascal} Pa.");

        Temperature Tsat = SaturationTemperature.Temperature(pTarget);

        // Ancillary-based rhoL/rhoV lose accuracy fast in the last ~1 K below Tc (see remarks
        // above) - reject rather than silently return an inaccurate two-phase state there.
        if (Tc - Tsat < CriticalMargin)
            throw new System.InvalidOperationException(
                $"UpdatePX: pressure too close to the critical point for a reliable fast ancillary-based " +
                $"two-phase solve (Tc-Tsat={(Tc - Tsat).Kelvin:G3} K). P={pTarget.Pascal} Pa. Use UpdatePT with a phase hint instead.");

        Molarity rhoL = BubbleDensity.Density(Tsat);
        Molarity rhoV = DewDensity.Density(Tsat);

        var sat = new SaturationSolver.SatResult(Tsat, pTarget, rhoL, rhoV);
        local.SetTwoPhase(sat, quality);
    }
}
