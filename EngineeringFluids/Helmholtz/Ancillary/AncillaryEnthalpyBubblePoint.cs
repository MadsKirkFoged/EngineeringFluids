using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz.Ancillary;
public static class AncillaryEnthalpyBubblePoint
{

    public static readonly List<double> A = [ -14548.431641457397,
                                              -425.3097032919506,
                                               6.040305720967621,
                                              -0.03667979513763329,
                                               0.00013121909101221487,
                                              -2.833690304020896e-07,
                                               3.4105973158767573e-10,
                                              -1.7632569242236018e-13];


    public static readonly List<double> B = [1, -0.002449612160719887];


    public static readonly double maxAbsError = 175.02638309168015;
    public static readonly double Tmax = 405.3;
    public static readonly double Tmin = 195.495;

    public static readonly double hs_anchor = 25012.280343289814;

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

        return numerator / denominator + hs_anchor;
    }

}
