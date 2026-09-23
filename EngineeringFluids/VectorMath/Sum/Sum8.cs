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


   
}
