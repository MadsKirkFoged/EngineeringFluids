using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class ResidualHelmholtzGaoB
{
    private static readonly double[] b = { 1.244, 0.6826 };
    private static readonly double[] beta = { 0.3696, 0.2962 };
    private static readonly double[] d = { 1, 1 };
    private static readonly double[] epsilon = { 0.4478, 0.44689 };
    private static readonly double[] eta = { -2.8452, -2.8342 };
    private static readonly double[] gamma = { 1.108, 1.313 };
    private static readonly double[] n = { -1.6909858, 0.93739074 };
    private static readonly double[] t = { 4.3315, 4.015 };

    public static double alphaR(double delta, double tau)
    {
        // Calculate alphaR
        double alphaRb = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {

            //Ftau = pow(tau, t) * exp(1.0 / (b + beta * pow(-gamma + tau, 2)));
            double Ftau = Math.Pow(tau, t[i]) * Math.Exp(1.0 / (b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2)));

            //Fdelta = pow(delta, d) * exp(eta * pow(delta - epsilon, 2));
            double Fdelta = Math.Pow(delta, d[i]) * Math.Exp(eta[i] * Math.Pow(delta - epsilon[i], 2));
            double term = n[i] * Ftau * Fdelta;
            alphaRb += term;
        }
        return alphaRb;
    }


    public static double alphaR_dDelta(double delta, double tau)
    {
        double dalphaRb_dDelta = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {

            double Ftau = Math.Pow(tau, t[i]) * Math.Exp(1.0f / (b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2)));

            double deltaTerm = d[i] * Math.Pow(delta, d[i] - 1);
            double expTerm = Math.Exp(eta[i] * Math.Pow(delta - epsilon[i], 2));
            double expTermDerivative = 2 * eta[i] * (delta - epsilon[i]) * expTerm;

            // Combined term for derivative:

            double term = n[i] * Ftau * (deltaTerm * expTerm + Math.Pow(delta, d[i]) * expTermDerivative);

            dalphaRb_dDelta += term;
        }

        return dalphaRb_dDelta;
    }
       


    //From coolprop
    public static double alphaR_dTau(double delta, double tau)
    {

        // Calculate alphaR
        double alphaRb = 0.0;
        for (int i = 0; i < n.Length; i++)
        {
            
            double Fdelta = Math.Pow(delta, d[i]) * Math.Exp(eta[i] * Math.Pow(delta - epsilon[i], 2));


            double taudFtaudtau = (2 * beta[i] * Math.Pow(tau, t[i] + 1) * (gamma[i] - tau) + t[i] * Math.Pow(tau, t[i]) * Math.Pow(b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2), 2))
                       * Math.Exp(1.0 / (b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2))) / Math.Pow(b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2), 2);

            double term = n[i] * Fdelta * taudFtaudtau / tau;
            

            alphaRb += term;

        }
        return alphaRb;
    }

    public static double alphaR_dDelta2(double delta, double tau)
    {
        double d2 = 0.0;

        for (int i = 0; i < n.Length; i++)
        {
            // Ftau = tau^t * exp( 1/(b + beta*(gamma - tau)^2) )
            double Ftau = Math.Pow(tau, t[i]) *
                          Math.Exp(1.0 / (b[i] + beta[i] * Math.Pow(gamma[i] - tau, 2)));

            // Fdelta = delta^d * exp(eta*(delta - epsilon)^2)
            // Need second derivative of Fdelta
            int di = (int)d[i];

            double g = Math.Pow(delta, di);
            double gp = (di == 0) ? 0.0 : di * Math.Pow(delta, di - 1);
            double gpp = (di <= 1) ? 0.0 : di * (di - 1) * Math.Pow(delta, di - 2);

            // u(delta) = eta*(delta-epsilon)^2
            // u' = 2*eta*(delta-epsilon)
            // u'' = 2*eta
            double u1 = 2.0 * eta[i] * (delta - epsilon[i]);
            double u2 = 2.0 * eta[i];

            double E = Math.Exp(eta[i] * Math.Pow(delta - epsilon[i], 2));

            // Fdelta'' = E*[g'' + 2*g'*u' + g*(u'' + u'^2)]
            double Fdelta_pp = E * (gpp + 2.0 * gp * u1 + g * (u2 + u1 * u1));

            d2 += n[i] * Ftau * Fdelta_pp;
        }

        return d2;
    }



}
