using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Density + temperature: the EOS is explicit in (T, rho), so the only work is deciding whether the state is
    // inside the dome. Below Tc that uses the exact saturated densities at T (SolveAtT), like CoolProp: below the
    // vapor density it is gas, above the liquid density liquid, in between two-phase with the quality from the
    // lever rule on specific volume. A density within DomeEdgeBand of a saturated density counts as saturated, so a
    // state CoolProp puts exactly on the saturation line is found on it here too.
    //
    // Labels: CoolProp's plain single-phase rule (its DT flash labels are recalculated with it afterwards).
    public static void UpdateDT(this Ammonia local, Density density, Temperature temperature)
    {
        if (!double.IsFinite(density.KilogramPerCubicMeter) || density <= Density.Zero)
            throw new ArgumentOutOfRangeException(nameof(density), $"UpdateDT: density must be positive. D={density.KilogramPerCubicMeter} kg/m3.");
        if (!double.IsFinite(temperature.Kelvin) || temperature < local.TripleLiquid.Temperature)
            throw new ArgumentOutOfRangeException(nameof(temperature),
                $"UpdateDT: temperature below the triple point ({local.TripleLiquid.Temperature.Kelvin} K) is not supported. T={temperature.Kelvin} K.");

        local.ClearTwoPhase();
        SetDT(local, density, temperature);
    }

    // The DT flash itself, also the inner evaluator of the DH and DS flashes
    private static void SetDT(Ammonia local, Density density, Temperature T)
    {
        Temperature Tc = local.Critical.Temperature;
        Molarity rhoc = local.Critical.MolarDensity;
        Molarity rho = density / local.MolarMass;

        if (T < Tc && TrySaturationAtT(local, T) is SaturationSolver.SatResult sat
            && rho >= (1 - DomeEdgeBand) * sat.RhomolarV && rho <= (1 + DomeEdgeBand) * sat.RhomolarL)
        {
            double q = (double)((1.0 / rho - 1.0 / sat.RhomolarL) / (1.0 / sat.RhomolarV - 1.0 / sat.RhomolarL));
            local.SetTwoPhase(sat, Math.Clamp(q, 0.0, 1.0));
            return;
        }

        local.ClearTwoPhase();
        local.Temperature = T;
        local.Density = density;
        if (T == Tc && rho == rhoc)
            local.SetPhase(Phases.CriticalPoint);
    }

    // Exact saturation at T, including the triple point. SolveAtT converges to within 6 microK of Tc (tested); if it
    // fails closer than that, the dome is too narrow to resolve and the state is treated as single-phase (null).
    private static SaturationSolver.SatResult? TrySaturationAtT(Ammonia local, Temperature T)
    {
        try
        {
            return local.SaturationAtT(T);
        }
        catch (InvalidOperationException)
        {
            return null;
        }
    }
}
