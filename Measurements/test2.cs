using System;
using System.Linq;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Diagnosers;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Running;
using EngineeringFluids.Helmholtz;

[MemoryDiagnoser]
[DisassemblyDiagnoser(maxDepth: 3)]
[ShortRunJob]
// Uncomment if you want to see inlining decisions
// [InliningDiagnoser]
[SimpleJob(RuntimeMoniker.Net90, baseline: true)]
public class PressureBenchmarks
{
    private float[] _randomTemps = default!;
    private float[] _rampTemps = default!;
    private int _idx;

    // Controls how many distinct inputs we cycle through.
    // More inputs reduces caching/branch prediction artifacts.
    [Params(256, 4096)]
    public int InputCount;

    // Used by the tight-loop benchmark
    [Params(64, 1024)]
    public int LoopCount;

    [GlobalSetup]
    public void Setup()
    {
        var rng = new Random(42);
        _randomTemps = new float[InputCount];
        _rampTemps = new float[InputCount];

        float tMin = 195.495f;
        float tMax = 405.56f;

        for (int i = 0; i < InputCount; i++)
        {
            // Uniform random in [tMin, tMax]
            _randomTemps[i] = tMin + (float)rng.NextDouble() * (tMax - tMin);

            // Smooth ramp (monotonic) across the domain
            _rampTemps[i] = tMin + (tMax - tMin) * (i / (float)(InputCount - 1));
        }

        _idx = 0;
    }

    [Benchmark(Baseline = true)]
    public float OldCorrelation_RandomScalar()
    {
        float t = NextRandomTemp();
        return Saturation.CalculateSaturationPressureFloat(t);
    }

    [Benchmark]
    public float FastPressurePoly_RandomScalar()
    {
        float t = NextRandomTemp();
        return BubbleEnthalpyFast.BubbleEnthalpy(t);
    }

    //[Benchmark]
    //public float TwoPartQuadratic_RandomScalar()
    //{
    //    float t = NextRandomTemp();
    //    return SaturationPressureFast.Pressure(t);
    //}

    //[Benchmark]
    //public float OldCorrelation_RampScalar()
    //{
    //    float t = NextRampTemp();
    //    return Saturation.CalculateSaturationPressureFloat(t);
    //}

    //[Benchmark]
    //public float FastPressurePoly_RampScalar()
    //{
    //    float t = NextRampTemp();
    //    return FastPressurePoly.Pressure(t);
    //}

    //[Benchmark]
    //public float TwoPartQuadratic_RampScalar()
    //{
    //    float t = NextRampTemp();
    //    return SaturationPressureFast.Pressure(t);
    //}

    // Tight loop throughput benchmark: simulates many calls in a tight loop
    // (more sensitive to branch predictability / cache / instruction mix).
    //[Benchmark]
    //public float TwoPartQuadratic_LoopRandom()
    //{
    //    float acc = 0;
    //    for (int i = 0; i < LoopCount; i++)
    //        acc += SaturationPressureFast.Pressure(NextRandomTemp());
    //    return acc;
    //}

    //[Benchmark]
    //public float FastPressurePoly_LoopRandom()
    //{
    //    float acc = 0;
    //    for (int i = 0; i < LoopCount; i++)
    //        acc += FastPressurePoly.Pressure(NextRandomTemp());
    //    return acc;
    //}

    private float NextRandomTemp()
    {
        // Bitmask is faster than modulo when InputCount is a power of two.
        // If you use non-pow2 counts, change to: _idx = (_idx + 1) % InputCount;
        _idx++;
        return _randomTemps[_idx & (InputCount - 1)];
    }

    private float NextRampTemp()
    {
        _idx++;
        return _rampTemps[_idx & (InputCount - 1)];
    }
}

//public static class Program
//{
//    public static void Main(string[] args)
//        => BenchmarkRunner.Run<PressureBenchmarks>();
//}
