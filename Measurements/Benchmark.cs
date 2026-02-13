using BenchmarkDotNet.Attributes;
using EngineeringFluids;
using System.Numerics;
using SharpFluids;
using EngineeringFluids.Fluids;
using EngineeringUnits;

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
    private static readonly Ammonia EngineeringFluid = new Ammonia();

    private static readonly Pressure pressure = Pressure.FromBar(21);
    private static readonly Temperature temperature = Temperature.FromDegreeCelsius(100);



    // 1. speed upgrade
    //    | Method              | Mean      | Error     | StdDev    | Gen0    | Allocated |
    //|-------------------- |----------:|----------:|----------:|--------:|----------:|
    //| SharpFluidupdatePT  |  22.78 us |  10.28 us |  0.563 us |  0.7629 |   3.26 KB |
    //| EngineeringUpdatePT | 473.09 us | 185.31 us | 10.157 us | 70.3125 | 299.65 KB |

    // 2. speed upgrade (It make it slower...)
    //    | Method              | Mean        | Error     | StdDev    | Gen0     | Allocated |
    //|-------------------- |------------:|----------:|----------:|---------:|----------:|
    //| SharpFluidupdatePT  |    21.38 us |  32.44 us |  1.778 us |   0.7629 |   3.26 KB |
    //| EngineeringUpdatePT | 1,227.37 us | 448.49 us | 24.583 us | 289.0625 |   1225 KB |


    [Benchmark]
    public void SharpFluidupdatePT()
    {
        SharpFluid.UpdatePT(pressure, temperature);
    }


    [Benchmark]
    public void EngineeringUpdatePT()
    {
        EngineeringFluid.UpdatePT(pressure, temperature);
    }



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