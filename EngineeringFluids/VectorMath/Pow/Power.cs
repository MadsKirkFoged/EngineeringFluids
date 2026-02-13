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
public static class Power
{

    public static float Pow(float x, int y)
    {
        var result = 1f;

        for (int i = 0; i < y; i++)
        {
            result *= x;
        }
        return result;
    }

}
