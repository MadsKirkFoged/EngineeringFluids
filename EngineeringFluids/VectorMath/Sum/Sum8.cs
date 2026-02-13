using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Reflection;


//using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;



namespace EngineeringMath.VectorMath.Pow;
public static class Sum8
{

    //Fastest but almost the same speed as just adding them all together
    public static void Sum(Span<float> result)
    {
        Vector256<float> v1 = new Vector<float>(result).AsVector256();

        result[0] = Vector256.Sum(v1);
    }


    //public static void SumGPT2(Span<float> span)
    //{
    //    span[0] = span[0] + span[1] + span[2] + span[3] + span[4] + span[5] + span[6] + span[7];
    //}


    //public static void SumGPT(Span<float> result)
    //{
    //    var vector = new Vector<float>(result);

    //    result[0] = Vector.Dot(vector, Vector<float>.One);
    //}


    //public static void SumGPT3(Span<float> span)
    //{
    //    float sum = 0;
    //    for (int i = 0; i < span.Length; i++)
    //    {
    //        sum += span[i];
    //    }

    //    span[0] = sum;
    //}

}
//private static float[] Add1 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
//private static float[] Add2 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
//private static float[] Add3 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
//private static float[] Add4 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];
//private static float[] Add5 = [0.000567f, 12345678.9f, 34567890.1f, 45678901.2f, 0.001234f, 54353.343467f, 12345678.9f, 4f];


//[Benchmark]
//public static void sum() => Sum8.Sum(Add1);

//[Benchmark]
//public static void SumGPT() => Sum8.SumGPT(Add2);

//[Benchmark]
//public static void SumGPT2() => Sum8.SumGPT2(Add3);

//[Benchmark]
//public static void SumGPT3() => Sum8.SumGPT3(Add4);