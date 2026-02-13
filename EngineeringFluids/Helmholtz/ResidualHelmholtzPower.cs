using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class ResidualHelmholtzPower
{

    
    private static readonly double[] n = { 0.006132232, 1.7395866, -2.2261792, -0.30127553, 0.08967023, -0.076387037, -0.84063963, -0.27026327 };
    private static readonly int[] d = { 4, 1, 1, 2, 3, 3, 2, 3 };
    private static readonly double[] t_R = { 1.0, 0.382, 1.0, 1.0, 0.677, 2.915, 3.51, 1.063 };
    private static readonly int[] l = { 0, 0, 0, 0, 0, 2, 2, 1 };



    public static double alphaR(double delta, double tau)
    {
        double alphaRP = 0.0;


        for (int i = 0; i < n.Length; i++)
        {
            double term = n[i] * EngineeringMath.Pow(delta, d[i]) * EngineeringMath.PowWithInts(tau, t_R[i]);

            if (l[i] != 0)
            {
                term *= Math.Exp(-EngineeringMath.Pow(delta, l[i]));
            }
            alphaRP += term;
        }

        return alphaRP;
    }


    public static double alphaR_dDelta(double delta, double tau)
    {
        double dalphaRP_dDelta = 0.0;

        for (int i = 0; i < n.Length; i++)
        {
            // Base term without the exponential component
            double baseTerm = n[i] * EngineeringMath.PowWithInts(tau, t_R[i]);
            // Derivative of the power term
            double powerTerm = d[i] * EngineeringMath.Pow(delta, d[i] - 1);
            // Initial term calculation
            double term = baseTerm * powerTerm;

            if (l[i] != 0)
            {
                double expTerm = Math.Exp(-EngineeringMath.Pow(delta, l[i]));
                double expDerivativeTerm = -l[i] * EngineeringMath.Pow(delta, l[i] - 1) * expTerm;

                term = baseTerm * (powerTerm * expTerm + EngineeringMath.Pow(delta, d[i]) * expDerivativeTerm);
            }

            dalphaRP_dDelta += term;
        }

        return dalphaRP_dDelta;
    }


    public static double alphaR_dTau(double delta, double tau)
    {
        double dalphaRP_dTau = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {
            double term = n[i] * Math.Pow(delta, d[i]) * Math.Pow(tau, t_R[i]);

            if (l[i] != 0)
            {
                term *= Math.Exp(-Math.Pow(delta, l[i]));
            }

            // Differentiate term with respect to tau
            double dterm_dTau = term * t_R[i] / tau;

            dalphaRP_dTau += dterm_dTau;
        }

        return dalphaRP_dTau;
    }

    public static double alphaR2_dTau(double delta, double tau)
    {
        double alphaRP = 0.0;


        for (int i = 0; i < n.Length; i++)
        {
            double term = n[i] * EngineeringMath.Pow(delta, d[i]) * t_R[i] * EngineeringMath.PowWithInts(tau, t_R[i] - 1);

            if (l[i] != 0)
            {
                term *= Math.Exp(-EngineeringMath.Pow(delta, l[i]));
            }
            alphaRP += term;
        }

        return alphaRP;
    }

    public static double alphaR_dDelta2(double delta, double tau)
    {
        double d2 = 0.0;

        for (int i = 0; i < n.Length; i++)
        {
            double tauPart = EngineeringMath.PowWithInts(tau, t_R[i]);
            double A = n[i] * tauPart;

            int di = d[i];
            int li = l[i];

            if (li == 0)
            {
                // term = A * delta^d
                // term'' = A * d*(d-1) * delta^(d-2)
                if (di >= 2)
                {
                    d2 += A * di * (di - 1) * EngineeringMath.Pow(delta, di - 2);
                }
                continue;
            }

            // term = A * g * E
            // g = delta^d
            // E = exp(-delta^l)

            double g = EngineeringMath.Pow(delta, di);
            double E = Math.Exp(-EngineeringMath.Pow(delta, li));

            double gp = (di == 0) ? 0.0 : di * EngineeringMath.Pow(delta, di - 1);
            double gpp = (di <= 1) ? 0.0 : di * (di - 1) * EngineeringMath.Pow(delta, di - 2);

            // Useful powers (all non-negative because li >= 1 here)
            double delta_lm1 = EngineeringMath.Pow(delta, li - 1);           // delta^(l-1)
            double delta_2l_2 = EngineeringMath.Pow(delta, 2 * li - 2);      // delta^(2l-2)
            double delta_lm2 = (li >= 2) ? EngineeringMath.Pow(delta, li - 2) : 0.0; // delta^(l-2) only if l>=2

            // Derived formula:
            // term'' = A*E*[ g'' - 2*l*delta^(l-1)*g' + g*(l^2*delta^(2l-2) - l*(l-1)*delta^(l-2)) ]
            double bracket =
                gpp
                - 2.0 * li * delta_lm1 * gp
                + g * (li * li * delta_2l_2 - li * (li - 1) * delta_lm2);

            d2 += A * E * bracket;
        }

        return d2;
    }




}
