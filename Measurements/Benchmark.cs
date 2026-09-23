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




}