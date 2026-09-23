using BenchmarkDotNet.Attributes;
using EngineeringFluids;
using System.Numerics;
using System.Numerics.Tensors;
using SharpFluids;
using EngineeringFluids.Fluids;
using EngineeringUnits;
using EngineeringFluids.Helmholtz;

namespace Benchmark;

[ShortRunJob]
[MemoryDiagnoser]
public class Benchy
{
    private static readonly double T = 400.0; // Temperature in K
    private static readonly double rho = 9.0; // Density in kg/m^3


    private static readonly double Tc = 405.56; // Critical temperature in K
    private static readonly double rhoc = 13696.0; // Critical density in mol/m^3
    private static readonly double M = 0.01703056; // Molar mass in kg/mol
    private static readonly double R = 8.3144598; // Gas constant in J/(mol*K)

    private static readonly Pressure pressure = Pressure.FromBar(21);
    private static readonly Temperature temperature = Temperature.FromDegreeCelsius(100);

    private static readonly double pressured = 2100000; // Pressure.FromBar(21).SI;
    private static readonly double temperatured = 373.15; // Temperature.FromDegreeCelsius(100).SI;

    private static readonly Fluid SharpFluid = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluid = new AmmoniaDouble();
    private static readonly AbstractState CoolPropState = AbstractState.factory("HEOS", "Ammonia");

    [Benchmark]
    public void SharpFluidupdatePT()
    {
        SharpFluid.UpdatePT(pressure, temperature);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePT()
    {
        CoolPropState.update(input_pairs.PT_INPUTS, pressured, temperatured);
    }

    [Benchmark]
    public void EngineeringUpdatePT()
    {
        EngineeringFluid.UpdatePT(pressured, temperatured);
    }

    // ---------------------------------------------------------------------
    // Region-specific cases for the Fast (AmmoniaDouble) path only.
    //
    // UpdatePT's Newton solver needs a different number of iterations depending
    // on where (T,P) falls relative to the saturation dome and the critical
    // point (measured via a full-domain sweep while tuning the initial-guess
    // logic): ~2-3 iterations for superheated gas far from the dome, ~4-5 near
    // it or in the liquid region, and up to an order of magnitude more right
    // near the critical point. The single point above (superheated gas, far
    // from the dome) only exercises the cheapest case, so a regression that
    // only shows up in liquid/supercritical/near-critical states would be
    // invisible there. These cases give each region its own tracked number.
    // ---------------------------------------------------------------------

    private static readonly AmmoniaDouble GasNearDomeFluid = new AmmoniaDouble();
    private static readonly double gasNearDomePressure = 5_630_000; // ~0.9x Psat(373.15K): gas, close to the dew line
    private static readonly double gasNearDomeTemperature = 373.15;

    private static readonly AmmoniaDouble LiquidFluid = new AmmoniaDouble();
    private static readonly double liquidPressure = 2_100_000;
    private static readonly double liquidTemperature = 280.0; // subcooled liquid

    private static readonly AmmoniaDouble SupercriticalTypicalFluid = new AmmoniaDouble();
    private static readonly double supercriticalTypicalPressure = 15_000_000; // ~1.3x Pc, comfortably away from Tc
    private static readonly double supercriticalTypicalTemperature = 500.0;

    private static readonly AmmoniaDouble SupercriticalNearCriticalFluid = new AmmoniaDouble();
    private static readonly double supercriticalNearCriticalPressure = 14_158_915.68768298; // worst case found while sweeping the supercritical region
    private static readonly double supercriticalNearCriticalTemperature = 406.56; // 1 K above Tc=405.56K

    [Benchmark]
    public void Gas_NearDome_Fast()
    {
        GasNearDomeFluid.UpdatePT(gasNearDomePressure, gasNearDomeTemperature);
    }

    [Benchmark]
    public void Liquid_Fast()
    {
        LiquidFluid.UpdatePT(liquidPressure, liquidTemperature);
    }

    [Benchmark]
    public void Supercritical_Typical_Fast()
    {
        SupercriticalTypicalFluid.UpdatePT(supercriticalTypicalPressure, supercriticalTypicalTemperature);
    }

    [Benchmark]
    public void Supercritical_NearCritical_Fast()
    {
        SupercriticalNearCriticalFluid.UpdatePT(supercriticalNearCriticalPressure, supercriticalNearCriticalTemperature);
    }
}