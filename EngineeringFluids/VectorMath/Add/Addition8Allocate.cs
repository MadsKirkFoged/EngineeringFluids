using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace EngineeringMath.VectorMath.Add;

public static class Addition8Allocate
{

    //Test results:

    //Reuse -->
    // Add x + x   -->  0.4900 ns
    // AddVecFixed -->  1.4954 ns

    //No Reuse -->
    // AddFixed8      -->  8.476 ns

    public static float[] AddFixed8(ReadOnlySpan<float> Input1, ReadOnlySpan<float> Input2)
    {
        //Vector<float>.Count needs to be 8
        
        float[] Result = new float[8];

        Vector<float> v1 = new(Input1);
        Vector<float> v2 = new(Input2);

        (v1 + v2).CopyTo(Result);

        return Result;

    }


    // | Method        | Mean     | Error     | StdDev    | Gen0   | Allocated |
    // |-------------- |---------:|----------:|----------:|-------:|----------:|
    // | BenchMarkAdd1 | 8.476 ns | 0.2125 ns | 0.2087 ns | 0.0130 |      56 B |

    //private static float[] Add1 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
    //private static float[] Add2 = [101234567.8f, 67890123.4f, 0.001012f, 11f, 8345.435f, 45678901.2f, 0.001234f, 89012345.6f];

    //[Benchmark]
    //public static float[] BenchMarkAdd1() => Addition8Allocate.AddFixed8(Add1, Add2);

}
