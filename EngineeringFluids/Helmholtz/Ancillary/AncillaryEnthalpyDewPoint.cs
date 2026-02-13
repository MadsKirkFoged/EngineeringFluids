using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz.Ancillary;
public static class AncillaryEnthalpyDewPoint
{

    public static readonly List<double> A = [ -39656.441243327485,
                                              1732.1248720835576,
                                              -20.184954966702787,
                                              0.12336506690457887,
                                              -0.00044866664143201863,
                                              9.727626573853517e-07,
                                              -1.1638181243435474e-09,
                                              5.938354498744747e-13 ];


    public static readonly List<double> B = [1, -0.00244523421040168];


    public static readonly double maxAbsError = 360.2372587421132;
    public static readonly double Tmax = 405.3;
    public static readonly double Tmin = 195.495;

    public static readonly double hs_anchor = 25012.280343289814;


    //Output J/mol

    public static double MolarEnthalpy(Temperature T)
    {
        double Temp = T.Kelvin;

        double numerator = 0;
        for (int i = 0; i < A.Count; i++)
        {
            numerator += A[i] * Math.Pow(Temp, i);
        }

        double denominator = 0;
        for (int i = 0; i < B.Count; i++)
        {
            denominator += B[i] * Math.Pow(Temp, i);
        }

        return (numerator / denominator);// + hs_anchor;
    }

}
