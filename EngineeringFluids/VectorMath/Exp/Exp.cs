using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringMath.VectorMath.Exp;
public static class ExpClass
{


    public static void Exp8(Span<float> result)
    {
        Vector<float> v1 = new(result);

        Vector256<float> v2 = v1.AsVector256();
        Vector256<float> test = Vector256.Exp(v2);

        test.CopyTo(result);
    }

   
}
