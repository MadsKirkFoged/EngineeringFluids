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
}