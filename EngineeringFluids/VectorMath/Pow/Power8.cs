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
public static class Power8
{

    //  | Method   | Mean      | Error     | StdDev    | Allocated |
    //  |--------- |----------:|----------:|----------:|----------:|






    //  | PowInt1  |  2.924 ns | 0.0851 ns | 0.0664 ns |         - |
    public static void Pow(Span<float> result, int power)
    {
        Vector<float> v1 = new(result);
        Vector<float> start = Vector<float>.One;

        for (int i = 0; i < power; i++)
        {
            start *= v1;
        }

        start.CopyTo(result);

    }

    public static void Pow(Span<float> result, float power)
    {
        Vector256<float> tau = new Vector<float>(result).AsVector256();

        var v1 = Vector256.Log(tau);
        var v2 = Vector256.Multiply(v1, power);
        var v3 = Vector256.Exp(v2);

        v3.CopyTo(result);

    }

    //  | PowFloat | 39.703 ns | 0.1254 ns | 0.1112 ns |         - |
    public static void Pow(Span<float> result, Span<float> power)
    {
        Vector256<float> tau = new Vector<float>(result).AsVector256();
        Vector256<float> p = new Vector<float>(power).AsVector256();

        var v1 = Vector256.Log(tau);
        var v2 = Vector256.Multiply(v1, p);
        var v3 = Vector256.Exp(v2);

        v3.CopyTo(result);
    }


    //  | PowInt2  | 10.773 ns | 0.1186 ns | 0.1110 ns |         - |
    public static void Pow(Span<float> result, Span<int> power)
    {

        for (int i = 0; i < result.Length; i++)
        {
            if (power[i] is 0)
            {
                result[i] = 1;
                continue;
            }

            float local = result[i];

            for (int j = 1; j < power[i]; j++)
            {
                result[i] = result[i] * local;
            }

        }
    }


    public static Span<float> Pow(float input, Span<int> power)
    {
        Span<float> result = new float[8];
        result.Fill(input);

        Pow(result, power);

        return result;
    }



}
