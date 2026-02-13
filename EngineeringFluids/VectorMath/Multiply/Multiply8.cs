using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace EngineeringMath.VectorMath.Add;

public static class Multiply8
{

    //Test results:

    //Reuse -->
    // Add x + x      -->  0.4900 ns
    // AddFixed8Reuse -->  1.975 ns

    //No Reuse -->
    // AddFixed8      -->  8.055 ns



    public static void MultiplyFixed8Reuse(Span<float> result, ReadOnlySpan<float> Fixed)
    {
        //Vector<float>.Count needs to be 8

        Vector<float> v1 = new(result);
        Vector<float> v2 = new(Fixed);

        (v1 * v2).CopyTo(result);        
    }

    public static void MultiplyFixed8Reuse(Span<float> result, ReadOnlySpan<float> Fixed1, ReadOnlySpan<float> Fixed2)
    {
        //Vector<float>.Count needs to be 8

        Vector<float> v1 = new(result);
        Vector<float> v2 = new(Fixed1);
        Vector<float> v3 = new(Fixed2);

        (v1 * v2 * v3).CopyTo(result);
    }

    public static void MultiplyFixed8Reuse(Span<float> result, ReadOnlySpan<float> Fixed1, ReadOnlySpan<float> Fixed2, ReadOnlySpan<float> Fixed3)
    {
        //Vector<float>.Count needs to be 8

        Vector<float> v1 = new(result);
        Vector<float> v2 = new(Fixed1);
        Vector<float> v3 = new(Fixed2);
        Vector<float> v4 = new(Fixed3);

        (v1 * v2 * v3 * v4).CopyTo(result);
    }


    // | Method        | Mean     | Error     | StdDev    | Allocated |
    // |-------------- |---------:|----------:|----------:|----------:|
    // | BenchMarkAdd1 | 1.975 ns | 0.0571 ns | 0.0506 ns |         - |



    //private static float[] Add1 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
    //private static float[] Add2 = [101234567.8f, 67890123.4f, 0.001012f, 11f, 8345.435f, 45678901.2f, 0.001234f, 89012345.6f];

    //[Benchmark]
    //public static void BenchMarkAdd1() => Addition8.AddFixed8Reuse(Add1, Add2);



}
