using BenchmarkDotNet.Attributes;
using EngineeringFluids;
using System.Numerics;
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

    // Calculate reduced temperature and density
    private static readonly double tau = Tc / T;
    private static readonly double delta = rho / (rhoc * M);

    private static readonly float tauF = (float)(Tc / T);
    private static readonly float deltaF = (float)(rho / (rhoc * M));

    private static readonly Fluid SharpFluid = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluid = new AmmoniaDouble();

    private static readonly Pressure pressure = Pressure.FromBar(21);
    private static readonly Temperature temperature = Temperature.FromDegreeCelsius(100);

    private static readonly double pressured = 2100000; // Pressure.FromBar(21).SI;
    private static readonly double temperatured = 373.15; // Temperature.FromDegreeCelsius(100).SI;



    // 1. speed upgrade
    //    | Method              | Mean      | Error     | StdDev    | Gen0    | Allocated |
    //|-------------------- |----------:|----------:|----------:|--------:|----------:|
    //| SharpFluidupdatePT  |  22.78 us |  10.28 us |  0.563 us |  0.7629 |   3.26 KB |
    //| EngineeringUpdatePT | 473.09 us | 185.31 us | 10.157 us | 70.3125 | 299.65 KB |

    // 2. speed upgrade (After converting to double)
    //| Method              | Mean     | Error     | StdDev   | Gen0   | Allocated |
    //|-------------------- |---------:|----------:|---------:|-------:|----------:|
    //| SharpFluidupdatePT  | 20.30 us |  3.776 us | 0.207 us | 0.7629 |   3.26 KB |
    //| EngineeringUpdatePT | 83.37 us | 19.059 us | 1.045 us | 1.4648 |   6.45 KB |

    // Added a cache
    //    | Method              | Mean     | Error     | StdDev   | Gen0   | Allocated |
    //|-------------------- |---------:|----------:|---------:|-------:|----------:|
    //| SharpFluidupdatePT  | 20.68 us | 10.163 us | 0.557 us | 0.7629 |   3.26 KB |
    //| EngineeringUpdatePT | 11.00 us |  2.389 us | 0.131 us | 1.5259 |   6.45 KB |


    // No more units
    //    | Method              | Mean      | Error     | StdDev   | Gen0   | Allocated |
    //|-------------------- |----------:|----------:|---------:|-------:|----------:|
    //| SharpFluidupdatePT  | 26.635 us | 19.343 us | 1.060 us | 0.7629 |   3.26 KB |
    //| EngineeringUpdatePT |  7.353 us | 32.568 us | 1.785 us | 0.9003 |   3.81 KB |


    //No unit and more cache
    //    | Method              | Mean        | Error       | StdDev    | Gen0   | Allocated |
    //|-------------------- |------------:|------------:|----------:|-------:|----------:|
    //| SharpFluidupdatePT  | 20,619.5 ns | 6,844.49 ns | 375.17 ns | 0.7629 |    3336 B |
    //| EngineeringUpdatePT |    990.1 ns |    63.84 ns |   3.50 ns |      - |         - |

//    | Method              | Mean     | Error    | StdDev   | Gen0   | Allocated |
//|-------------------- |---------:|---------:|---------:|-------:|----------:|
//| SharpFluidupdatePT  | 19.65 us | 4.600 us | 0.252 us | 0.7629 |    3336 B |
//| EngineeringUpdatePT | 42.14 us | 2.017 us | 0.111 us |      - |         - |

//    | Method              | Mean     | Error    | StdDev   | Gen0   | Allocated |
//|-------------------- |---------:|---------:|---------:|-------:|----------:|
//| SharpFluidupdatePT  | 19.57 us | 2.507 us | 0.137 us | 0.7629 |    3336 B |
//| EngineeringUpdatePT | 18.92 us | 0.916 us | 0.050 us |      - |         - |


    [Benchmark]
    public void SharpFluidupdatePT()
    {
        SharpFluid.UpdatePT(pressure, temperature);
    }


    [Benchmark]
    public void EngineeringUpdatePT()
    {
        EngineeringFluid.UpdatePT(pressured, temperatured);
    }





    //[Benchmark]
    //public double CalculateSaturationPressureDouble()
    //{
    //    return Saturation.CalculateSaturationPressureDouble(373.15d);
    //}

    //[Benchmark]
    //public float CalculateSaturationPressureFLoat()
    //{
    //    return Saturation.CalculateSaturationPressureFloat(373.15f);
    //}

    //[Benchmark]
    //public float FastPressurePolytest()
    //{
    //    return FastPressurePoly.Pressure(373.15f);
    //}

    //[Benchmark]
    //public float SaturationPressureFasttest()
    //{
    //    return SaturationPressureFast.Pressure(373.15f);
    //}

    //[Benchmark]
    //public double Alpha0() => IdealHelmholtzPlanckEinstein.Alpha0(delta, tau);
    //[Benchmark]
    //public double Alpha0Fast() => IdealHelmholtzPlanckEinsteinFast.Alpha0(delta, tau);


    //[Benchmark]
    //public double Alpha0_dTau() => IdealHelmholtzPlanckEinstein.Alpha0_dTau(delta, tau);
    //[Benchmark]
    //public double Alpha0_dTauFast() => IdealHelmholtzPlanckEinsteinFast.Alpha0_dTau(delta, tau);






    //| Method     | Mean        | Error     | StdDev    | Allocated |
    //|----------- |------------:|----------:|----------:|----------:|
    //| alpha0     |    52.53 ns |  0.527 ns |  0.440 ns |         - |
    //| alpha0Fast |    39.05 ns |  0.230 ns |  0.192 ns |         - |
    //| alphaR     | 1,567.16 ns | 15.581 ns | 14.574 ns |         - |



    //[Benchmark]
    //public static float alpha0() => AmmoniaOld.alpha0(deltaF, tauF);
    //[Benchmark]
    //public static float alphaR() => AmmoniaOld.alphaR(deltaF, tauF);


    //[Benchmark]
    //public static float alpha0Fast() => AmmoniaFast.alpha0(deltaF, tauF);
    //[Benchmark]
    //public static float alphaRFast() => AmmoniaFast.alphaR(deltaF, tauF);



    //[Benchmark]
    //public static double alpha0Float() => Ammonia.alpha0Float(deltaF, tauF);

    //[Benchmark]
    //public static double alpha0() => Ammonia.alpha0(delta, tau);

    //[Benchmark]
    //public static double alphaR() => Ammonia.alphaR(delta, tau);

    //| Method                          | Mean     | Error   | StdDev  | Allocated |
    //|-------------------------------- |---------:|--------:|--------:|----------:|
    //| ResidualHelmholtzPower          | 414.7 ns | 3.06 ns | 2.39 ns |         - |
    //| ResidualHelmholtzPowerOptimized | 196.7 ns | 3.20 ns | 2.99 ns |         - |

    //[Benchmark]
    //public static double ResidualHelmholtzPower() => Ammonia.ResidualHelmholtzPower(delta, tau);

    //[Benchmark]
    //public static double ResidualHelmholtzPowerOptimized() => Ammonia.ResidualHelmholtzPowerOptimized(delta, tau);

    //[Benchmark]
    //public static double ResidualHelmholtzPowerOptimizedF() => Ammonia.ResidualHelmholtzPowerOptimized(deltaF, tauF);
    //[Benchmark]
    //public static double ResidualHelmholtzPowerVecOptimized() => Ammonia.ResidualHelmholtzPowerVecOptimized(deltaF, tauF);


    //[Benchmark]
    //public static double ResidualHelmholtzGaussian() => Ammonia.ResidualHelmholtzGaussian(delta, tau);


    //[Benchmark]
    //public static double ResidualHelmholtzGaoB() => Ammonia.ResidualHelmholtzGaoB(delta, tau);


    private static readonly float[] left = Enumerable.Range(0, 1001).Select(i => (float)i).ToArray();
    private static readonly float[] right = Enumerable.Range(0, 1001).Select(i => (float)i).ToArray();



    //[Benchmark]
    //public static float[] SumVectortest() => SumVector(left, right);


    //[Benchmark]
    //public static float[] Sumtest() => Sum(left, right);

    //[Benchmark]
    //public static float[] SumNoChecktest() => SumNoCheck(left, right);

    //[Benchmark]
    //public static Span<float> SumSpanTest() => SumSpan(left, right);


    public static float[] SumVector(float[] left, float[] right)
    {
        if (left is null)
        {
            throw new ArgumentNullException(nameof(left));
        }

        if (right is null)
        {
            throw new ArgumentNullException(nameof(right));
        }

        if (left.Length != right.Length)
        {
            throw new ArgumentException($"{nameof(left)} and {nameof(right)} are not the same length");
        }

        int length = left.Length;
        float[] result = new float[length];

        // Get the number of elements that can't be processed in the vector
        // NOTE: Vector<T>.Count is a JIT time constant and will get optimized accordingly
        int remaining = length % Vector<float>.Count;

        for (int i = 0; i < length - remaining; i += Vector<float>.Count)
        {
            var v1 = new Vector<float>(left, i);
            var v2 = new Vector<float>(right, i);
            (v1 + v2).CopyTo(result, i);
        }

        for (int i = length - remaining; i < length; i++)
        {
            result[i] = left[i] * right[i];
        }

        return result;
    }


    public static float[] Sum(float[] left, float[] right)
    {

        if (left is null)
        {
            throw new ArgumentNullException(nameof(left));
        }

        if (right is null)
        {
            throw new ArgumentNullException(nameof(right));
        }

        if (left.Length != right.Length)
        {
            throw new ArgumentException($"{nameof(left)} and {nameof(right)} are not the same length");
        }

        int length = left.Length;
        float[] result = new float[length];


        for (int i = 0; i < length; i++)
        {
            result[i] = left[i] * right[i];
        }

        return result;
    }

    public static float[] SumNoCheck(float[] left, float[] right)
    {

        int length = left.Length;
        float[] result = new float[length];


        for (int i = 0; i < length; i++)
        {
            result[i] = left[i] * right[i];
        }

        return result;
    }

    public static Span<float> SumSpan(Span<float> left, Span<float> right)
    {
        int length = left.Length;
        Span<float> result = new float[length];


        for (int i = 0; i < length; i++)
        {
            result[i] = left[i] * right[i];
        }

        return result;
    }

}