using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;

public static partial class Update
{
    // Density + pressure. Below Pc the exact saturation state at P decides first: a density between the saturated
    // vapor and liquid densities is two-phase at Tsat. Otherwise p(T) along the isochore is monotonic, and T is
    // found with Brent's method in a bracket that stays out of the dome:
    //   gas (rho < rhoV(P)):    Tsat(P) .. Tmax   (the isochore is superheated vapor above Tsat)
    //   liquid (rho > rhoL(P)): Tdome(rho) .. Tsat(P), Tdome being where the saturated-liquid density equals rho -
    //                           below it the isochore runs inside the dome, where the EOS is unstable
    //   P >= Pc:                above the dome (rho > rhoc: from Tdome, else from Tc - the critical isotherm is
    //                           below Pc for every rho < rhoc) up to Tmax
    // Labels: CoolProp's plain single-phase rule (the DP flash ends with recalculate_singlephase_phase).
    public static void UpdateDP(this Ammonia local, Density density, Pressure pressure)
    {
        if (!double.IsFinite(density.KilogramPerCubicMeter) || density <= Density.Zero)
            throw new ArgumentOutOfRangeException(nameof(density), $"UpdateDP: density must be positive. D={density.KilogramPerCubicMeter} kg/m3.");
        if (!double.IsFinite(pressure.Pascal) || pressure <= Pressure.Zero)
            throw new ArgumentOutOfRangeException(nameof(pressure), $"UpdateDP: pressure must be positive. P={pressure.Pascal} Pa.");

        local.ClearTwoPhase();

        Pressure Pc = local.Critical.Pressure;
        Pressure Ptriple = local.TripleLiquid.Pressure;
        Temperature Tc = local.Critical.Temperature;
        Temperature Ttriple = local.TripleLiquid.Temperature;
        Molarity rhoc = local.Critical.MolarDensity;
        Molarity rho = density / local.MolarMass;

        Temperature lo, hi;
        if (pressure >= Ptriple && pressure < Pc && TrySolveAtP(local, pressure) is SaturationSolver.SatResult sat)
        {
            if (rho >= (1 - DomeEdgeBand) * sat.RhomolarV && rho <= (1 + DomeEdgeBand) * sat.RhomolarL)
            {
                double q = (double)((1.0 / rho - 1.0 / sat.RhomolarL) / (1.0 / sat.RhomolarV - 1.0 / sat.RhomolarL));
                local.SetTwoPhase(sat, Math.Clamp(q, 0.0, 1.0));
                return;
            }

            // Liquid: the ancillary dome temperature can land at or above Tsat for a density just above rhoL, so
            // never start the bracket closer than 0.01 K below Tsat (a slightly stretched liquid there is still
            // below the target pressure).
            (lo, hi) = rho < sat.RhomolarV
                ? (sat.T, ExactFlashTmax)
                : (Temperature.Min(LiquidDomeTemperature(local, rho), sat.T - Temperature.FromKelvin(0.01)), sat.T);
        }
        else if (pressure >= Ptriple)
        {
            // At or above Pc - or so close below it that SolveAtP cannot resolve the dome: search above the dome
            lo = rho > rhoc ? LiquidDomeTemperature(local, rho) : Tc;
            hi = ExactFlashTmax;
        }
        else if (pressure < Ptriple && rho <= local.TripleVapor.MolarDensity)
        {
            // Below the triple-point pressure only vapor exists
            lo = Ttriple;
            hi = ExactFlashTmax;
        }
        else
        {
            throw new InvalidOperationException(
                $"UpdateDP: no fluid state at D={density.KilogramPerCubicMeter} kg/m3 and P={pressure.Pascal} Pa (below the triple-point pressure only vapor exists).");
        }

        double Residual(double TK)
        {
            local.Temperature = Temperature.FromKelvin(TK);
            local.Density = density;
            return (double)((local.Pressure - pressure) / pressure);
        }

        double flo = Residual(lo.Kelvin), fhi = Residual(hi.Kelvin);
        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"UpdateDP: could not bracket a solution in [{lo.Kelvin:G6} K, {hi.Kelvin:G6} K]. D={density.KilogramPerCubicMeter} kg/m3, P={pressure.Pascal} Pa.");

        double root = RootFinder.Brent(Residual, lo.Kelvin, hi.Kelvin, flo, fhi, xTol: 1e-12);
        if (!(Math.Abs(Residual(root)) <= 1e-9))
            throw new InvalidOperationException($"UpdateDP: did not converge. D={density.KilogramPerCubicMeter} kg/m3, P={pressure.Pascal} Pa.");
    }

    // Temperature at which the saturated-liquid density equals rho (the ancillary fit, inverted by bisection - it
    // only has to keep the bracket out of the dome, not be exact). Denser than the triple-point liquid: the triple point.
    private static Temperature LiquidDomeTemperature(Ammonia local, Molarity rho)
    {
        Temperature lo = local.TripleLiquid.Temperature;
        Temperature hi = local.Critical.Temperature;
        if (rho >= BubbleDensity.Density(lo))
            return lo;

        for (int i = 0; i < 60; i++)
        {
            Temperature mid = 0.5 * (lo + hi);
            if (BubbleDensity.Density(mid) > rho)
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    }
}
