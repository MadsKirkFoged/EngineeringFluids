using System;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace EngineeringMath.VectorMath.Add;

public static class Subtraction8Allocate
{

    //Test results:

    //Reuse -->
    //Same as add

    public static float[] SubFixed8(ReadOnlySpan<float> Input1, ReadOnlySpan<float> Input2)
    {
        //Vector<float>.Count needs to be 8
        
        float[] Result = new float[8];

        Vector<float> v1 = new(Input1);
        Vector<float> v2 = new(Input2);

        (v1 - v2).CopyTo(Result);

        return Result;

    } 

}
