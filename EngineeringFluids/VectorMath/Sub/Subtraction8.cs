using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace EngineeringMath.VectorMath.Add;

public static class Subtraction8
{

    //Test results:

    //Same as add



    public static void SubFixed8Reuse(Span<float> result, ReadOnlySpan<float> Fixed)
    {
        //Vector<float>.Count needs to be 8

        Vector<float> v1 = new(result);
        Vector<float> v2 = new(Fixed);

        (v1 - v2).CopyTo(result);
    }

}
