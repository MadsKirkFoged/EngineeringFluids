using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class IdealGasHelmholtzLead
{


    public static double Alpha0(double delta, double tau)
    {
        return Math.Log(delta) + -6.59406093943886 + 5.60101151987913 * tau;
    }

    public static double Alpha0_dTau(double delta, double tau)
    {
        // The derivative of the function with respect to tau
        return 5.60101151987913;
    }



}
