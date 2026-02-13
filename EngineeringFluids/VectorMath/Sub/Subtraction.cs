using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace EngineeringMath.VectorMath.Add;

public static class Subtraction
{

    //Test results:

    //Reuse
    // Add x + x    -->  0.4900 ns
    // AddVecFixed  -->  1.4954 ns

    //No Reuse 
    // AddFixed8    -->  8.055 ns



    //General purpose addition using SIMD
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
            result[i] = result[i] - Fixed[i];
        }

    }

}
