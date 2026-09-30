using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;

public static partial class Update
{
    // The ancillary-based two-phase flashes (UpdatePX/UpdateTX/UpdatePH/UpdatePS) are not
    // trusted within this distance below Tc - see the remarks on UpdatePX. (It was 1 K with the
    // float piecewise-quadratic ancillaries these replaced, which were 27% off near Tc.)
    private static readonly Temperature CriticalMargin = Temperature.FromKelvin(1e-3);

    // Pressure + vapor-quality flash: sets a two-phase saturated state
    // at pTarget with mass quality `quality` (0 = saturated liquid, 1 = saturated vapor).
    //
    // Unlike UpdatePXExact, this does NOT solve the phase equilibrium: it reads the saturation
    // temperature and the two branch densities straight from the ancillary tables
    // (SaturationTemperature / BubbleDensity / DewDensity). Downstream Enthalpy/Entropy/
    // InternalEnergy come from Ammonia's SetTwoPhase machinery, which evaluates the full EOS at
    // (Tsat, rhoL) and (Tsat, rhoV), so their accuracy follows from how well the tables track the
    // EOS dome. Measured against UpdatePXExact (p, rho, h, s, u): <= 7e-10 relative from the
    // triple point to 10 mK below Tc, 2e-8 between 10 mK and 1 mK. Closer than that rhoL and rhoV
    // are ill-conditioned even for the exact solver (a 5% perturbed start re-converges only to
    // ~1e-5), so that last millikelvin is rejected below rather than returning a state that is
    // ~2e-5 off.
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

        // Within the last millikelvin below Tc even exact rhoL/rhoV are ill-conditioned (see remarks
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

    // Exact (EOS-based) counterpart to UpdatePX: the saturation state comes from SolveAtP instead of the
    // ancillary fits, so it matches CoolProp to ~1e-10 and works right up to the critical point - including the
    // critical pressure itself, which (like CoolProp) gives the critical point for every quality.
    public static void UpdatePXExact(this Ammonia local, Pressure pTarget, double quality)
    {
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new System.ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        Pressure Pc = local.Critical.Pressure;
        Pressure Ptriple = local.TripleLiquid.Pressure;

        if (!double.IsFinite(pTarget.Pascal) || pTarget < Ptriple)
            throw new System.ArgumentOutOfRangeException(nameof(pTarget),
                $"UpdatePXExact: pressure below the triple-point pressure ({Ptriple.Pascal} Pa) has no liquid-vapor equilibrium. P={pTarget.Pascal} Pa.");

        if (pTarget > Pc)
            throw new System.InvalidOperationException(
                $"UpdatePXExact invalid above the critical pressure. P={pTarget.Pascal} Pa, Pc={Pc.Pascal} Pa.");

        if (pTarget == Pc)
        {
            local.SetTwoPhase(SaturationSolver.CriticalPoint(local), quality, EngineeringFluids.Helmholtz.Phase.Phases.CriticalPoint);
            return;
        }

        local.SetTwoPhase(local.SolveAtP(pTarget), quality);
    }
}
