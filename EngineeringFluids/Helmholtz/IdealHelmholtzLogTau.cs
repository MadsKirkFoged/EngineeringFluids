using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class IdealHelmholtzLogTau
{


    public static double Alpha0(double delta, double tau)
    {
        return 3 * Math.Log(tau);
    }


    public static double Alpha0_dTau(double delta, double tau)
    {
        // The derivative of the function with respect to tau
        return 3 / tau;
    }


}
