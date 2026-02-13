using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids;
public static class FastMath
{

    public static void Add(Span<float> result, ReadOnlySpan<float> Fixed)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            (v1 + v2).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] + Fixed[i];
        }

    }

    public static void Sub(Span<float> result, ReadOnlySpan<float> Fixed)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            (v1 - v2).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] + Fixed[i];
        }

    }

    public static void Multiply(Span<float> result, float Fixed1)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            (v1 * Fixed1).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] * Fixed1;
        }

    }

    public static void Multiply(Span<float> result, ReadOnlySpan<float> Fixed1)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed1.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            (v1 * v2).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] * Fixed1[i];
        }

    }

    public static void Multiply(Span<float> result, ReadOnlySpan<float> Fixed1, ReadOnlySpan<float> Fixed2)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed1.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            ReadOnlySpan<float> rightSlice2 = Fixed2.Slice(i, batchSize);
            Vector<float> v3 = new(rightSlice2);

            (v1 * v2 * v3).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] * Fixed1[i] * Fixed2[i];
        }

    }

    public static void Multiply(Span<float> result, ReadOnlySpan<float> Fixed1, ReadOnlySpan<float> Fixed2, ReadOnlySpan<float> Fixed3)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed1.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            ReadOnlySpan<float> rightSlice2 = Fixed2.Slice(i, batchSize);
            Vector<float> v3 = new(rightSlice2);

            ReadOnlySpan<float> rightSlice3 = Fixed3.Slice(i, batchSize);
            Vector<float> v4 = new(rightSlice3);

            (v1 * v2 * v3 * v4).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] * Fixed1[i] * Fixed2[i] * Fixed3[i];
        }

    }


    public static void Divide(Span<float> result, ReadOnlySpan<float> Fixed)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);

        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            ReadOnlySpan<float> rightSlice = Fixed.Slice(i, batchSize);
            Vector<float> v2 = new(rightSlice);

            (v1 / v2).CopyTo(leftSlice);
        }

        for (int i = upperLimit; i < length; i++)
        {
            result[i] = result[i] + Fixed[i];
        }

    }

    public static void Exp(Span<float> result)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);


        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            Vector256<float> v2 = v1.AsVector256();
            Vector256<float> test = Vector256.Exp(v2);

            test.CopyTo(leftSlice);

        }

    }

    public static void Log(Span<float> result)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);


        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            Vector256<float> v2 = v1.AsVector256();
            Vector256<float> test = Vector256.Log(v2);

            test.CopyTo(leftSlice);

        }

    }

    public static void Sum(Span<float> result)
    {

        int length = result.Length;
        var batchSize = Vector<float>.Count;
        var upperLimit = length & ~(batchSize - 1);


        for (int i = 0; i < upperLimit; i += batchSize)
        {
            Span<float> leftSlice = result.Slice(i, batchSize);
            Vector<float> v1 = new(leftSlice);

            Vector256<float> v2 = v1.AsVector256();
            Vector256<float> test = Vector256.CreateScalar((Vector256.Sum(v2)));

            test.CopyTo(leftSlice);

        }

    }


    public static void Pow(Span<float> result, int power)
    {
        if (Vector<float>.Count is 8)
        {

            Vector<float> v1 = new(result);
            Vector<float> start = Vector<float>.One;

            for (int i = 0; i < power; i++)
            {
                start *= v1;
            }

            start.CopyTo(result);
        }
    }


    public static void Pow(Span<float> result, Span<float> power)
    {
        if (Vector<float>.Count is 8)
        {

            Vector256<float> tau = new Vector<float>(result).AsVector256();
            Vector256<float> p = new Vector<float>(power).AsVector256();

            var v1 = Vector256.Log(tau);
            var v2 = Vector256.Multiply(v1, p);
            var v3 = Vector256.Exp(v2);

            v3.CopyTo(result);
        }
    }

    public static void Pow(Span<float> result, Span<int> power, int max)
    {
        if (Vector<float>.Count is 8)
        {

            Vector<float> v1 = new(result);
            Vector<int> v2 = new(power);
            Vector<float> start = new(result);



            for (int i = 1; i < max; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    if (v2.GetElement(j) is 0)
                    {
                        start = start.WithElement(j, 1f);
                    }

                    if (v2.GetElement(j) <= i)
                    {
                        v1 = v1.WithElement(j, 1f);
                    }
                }

                start *= v1;
            }

            start.CopyTo(result);
        }
    }

 

    public static void IdealHelmholtzPlanckEinstein(Span<float> Tau, ReadOnlySpan<float> n_PE, ReadOnlySpan<float> t_PE)
    {


        Vector256<float> tau = new Vector<float>(Tau).AsVector256();
        Vector256<float> n = new Vector<float>(n_PE).AsVector256();
        Vector256<float> t = new Vector<float>(t_PE).AsVector256();

        //2.224f * MathF.Log(1f + -1 * MathF.Exp(-4.0585856593352405f * tau));

        var f1 = Vector256.Multiply(t, tau);
        var f2 = Vector256.Exp(f1);
        var f3 = Vector256.Multiply(-1, f2);
        var f4 = Vector256.Add(Vector256<float>.One, f3);
        var f5 = Vector256.Log(f4);
        var f6 = Vector256.Multiply(n, f5);

        f6.CopyTo(Tau);

        for (int i = 3; i < Tau.Length; i++)
        {
            Tau[i] = 0;
        }

        var f7 = new Vector<float>(Tau).AsVector256();

        var f8 = Vector256.Sum(f7);
        var f9 = Vector256.CreateScalar(f8);

        f9.CopyTo(Tau);



    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void IdealHelmholtzPlanckEinstein2(Span<float> Tau, ReadOnlySpan<float> n_PE, ReadOnlySpan<float> t_PE)
    {


        Vector256<float> tau = new Vector<float>(Tau).AsVector256();
        Vector256<float> n = new Vector<float>(n_PE).AsVector256();
        Vector256<float> t = new Vector<float>(t_PE).AsVector256();

        //n * MathF.Log(1 + -1 * MathF.Exp(t * tau));

        var f1 = Vector256.Multiply(t, tau);
        var f2 = Vector256.Exp(f1);
        var f3 = Vector256.Multiply(-1, f2);
        var f4 = Vector256.Add(Vector256<float>.One, f3);
        var f5 = Vector256.Log(f4);
        var f6 = Vector256.Multiply(n, f5);

        f6.CopyTo(Tau);

        for (int i = 3; i < Tau.Length; i++)
        {
            Tau[i] = 0;
        }

        var f7 = new Vector<float>(Tau).AsVector256();

        var f8 = Vector256.Sum(f7);
        var f9 = Vector256.CreateScalar(f8);

        f9.CopyTo(Tau);

    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ResidualHelmholtzPower(Span<float> Tau, ReadOnlySpan<float> n_PE, ReadOnlySpan<float> t_PE)
    {


        Vector256<float> tau = new Vector<float>(Tau).AsVector256();
        Vector256<float> n = new Vector<float>(n_PE).AsVector256();
        Vector256<float> t = new Vector<float>(t_PE).AsVector256();

        //n * MathF.Log(1 + -1 * MathF.Exp(t * tau));

        var f1 = Vector256.Multiply(t, tau);
        var f2 = Vector256.Exp(f1);
        var f3 = Vector256.Multiply(-1, f2);
        var f4 = Vector256.Add(Vector256<float>.One, f3);
        var f5 = Vector256.Log(f4);
        var f6 = Vector256.Multiply(n, f5);

        f6.CopyTo(Tau);

        for (int i = 3; i < Tau.Length; i++)
        {
            Tau[i] = 0;
        }

        var f7 = new Vector<float>(Tau).AsVector256();

        var f8 = Vector256.Sum(f7);
        var f9 = Vector256.CreateScalar(f8);

        f9.CopyTo(Tau);
    }



}
