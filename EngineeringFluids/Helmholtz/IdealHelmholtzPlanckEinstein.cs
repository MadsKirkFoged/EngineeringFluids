using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class IdealHelmholtzPlanckEinstein
{

    private static readonly double[] n_PE = [2.224, 3.148, 0.9579];
    private static readonly double[] t_PE = [-4.0585856593352405, -9.776605187888352, -17.829667620080876];

    private static readonly double[] c_PE = [1, 1, 1];
    private static readonly double[] d_PE = [-1, -1, -1];



    public static double Alpha0(double delta, double tau)
    {
        double alpha0_PE = 0.0;
        for (int i = 0; i < n_PE.Length; i++)
        {
            alpha0_PE += n_PE[i] * Math.Log(c_PE[i] + d_PE[i] * Math.Exp(t_PE[i] * tau));
        }


        return alpha0_PE;
    }

    public static double Alpha0_dTau(double delta, double tau)
    {
        double dalpha0_PE_dTau = 0.0f;
        for (int i = 0; i < n_PE.Length; i++)
        {
            double numerator = n_PE[i] * d_PE[i] * t_PE[i] * Math.Exp(t_PE[i] * tau);
            double denominator = c_PE[i] + d_PE[i] * Math.Exp(t_PE[i] * tau);
            dalpha0_PE_dTau += numerator / denominator;
        }

        return dalpha0_PE_dTau;
    }

}
