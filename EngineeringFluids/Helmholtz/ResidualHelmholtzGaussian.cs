using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;
public static class ResidualHelmholtzGaussian
{

    private static readonly double[] beta = { 1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85 };
    private static readonly double[] d = { 1, 1, 1, 2, 2, 1, 3, 3, 1, 1 };
    private static readonly double[] epsilon = { -0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576 };
    private static readonly double[] eta = { 0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68 };
    private static readonly double[] gamma = { 1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277 };

    private static readonly double[] n = { 6.212578, -5.7844357, 2.4817542, -2.3739168, 0.01493697, -3.7749264, 0.0006254348, -1.7359e-05, -0.13462033, 0.07749072839 };
    private static readonly double[] t = { 0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13, 4.5, 1.0, 4.0 };




    public static double alphaR(double delta, double tau)
    {
        // Calculate alphaR
        double alphaRg = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            double deltaTerm = Math.Pow(delta, d[i]);
            double tauTerm = Math.Pow(tau, t[i]);
            double expTerm = Math.Exp(-eta[i] * Math.Pow(delta - epsilon[i], 2) - beta[i] * Math.Pow(tau - gamma[i], 2));
            double term = n[i] * deltaTerm * tauTerm * expTerm;
            alphaRg += term;
        }
        return alphaRg;
    }


    public static double alphaR_dDelta(double delta, double tau)
    {
        double dalphaRg_dDelta = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {

            double deltaTerm = d[i] * Math.Pow(delta, d[i] - 1);


            double tauTerm = Math.Pow(tau, t[i]);


            double expTerm = Math.Exp(-eta[i] * Math.Pow(delta - epsilon[i], 2) - beta[i] * Math.Pow(tau - gamma[i], 2));
            double expTermDerivative = -2 * eta[i] * (delta - epsilon[i]) * expTerm;

            // Combined term:
            double term = n[i] * (deltaTerm * tauTerm * expTerm + Math.Pow(delta, d[i]) * tauTerm * expTermDerivative);

            dalphaRg_dDelta += term;
        }

        return dalphaRg_dDelta;
    }

    public static double alphaR_dTau(double delta, double tau)
    {
        double dalphaRg_dTau = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            double deltaTerm = Math.Pow(delta, d[i]);
            double tauTerm = Math.Pow(tau, t[i]);
            double expTerm = Math.Exp(-eta[i] * Math.Pow(delta - epsilon[i], 2) - beta[i] * Math.Pow(tau - gamma[i], 2));
            double term = n[i] * deltaTerm * tauTerm * expTerm;

            // Differentiate term with respect to tau
            double dterm_dTau = term * (t[i] / tau - 2 * beta[i] * (tau - gamma[i]));

            dalphaRg_dTau += dterm_dTau;
        }
        return dalphaRg_dTau;
    }

    public static double alphaR2_dTau(double delta, double tau)
    {
        // Calculate alphaR
        double alphaRg = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            double deltaTerm = Math.Pow(delta, d[i]); //Same

            //double tauTerm = Math.Pow(tau, t_Rg[i]); //Gone

            double expTerm = Math.Exp(-eta[i] * Math.Pow(delta - epsilon[i], 2) - beta[i] * Math.Pow(tau - gamma[i], 2));

            double NewTerm = t[i] * Math.Pow(tau, t[i] -1) - 2 * beta[i] * Math.Pow(tau, t[i]) * (tau - gamma[i]);

            double term = n[i] * deltaTerm * NewTerm * expTerm;
            alphaRg += term;
        }
        return alphaRg;
    }

    public static double alphaR_dDelta2(double delta, double tau)
    {
        double d2 = 0.0;

        for (int i = 0; i < n.Length; i++)
        {
            int di = (int)d[i];

            double deltaPow = Math.Pow(delta, di);
            double tauPow = Math.Pow(tau, t[i]);

            double expAll = Math.Exp(
                -eta[i] * Math.Pow(delta - epsilon[i], 2)
                - beta[i] * Math.Pow(tau - gamma[i], 2)
            );

            double gp = (di == 0) ? 0.0 : di * Math.Pow(delta, di - 1);
            double gpp = (di <= 1) ? 0.0 : di * (di - 1) * Math.Pow(delta, di - 2);

            // u(delta) = -eta*(delta-eps)^2
            // u' = -2*eta*(delta-eps)
            // u'' = -2*eta
            double u1 = -2.0 * eta[i] * (delta - epsilon[i]);
            double u2 = -2.0 * eta[i];

            // d2(term) = n*tau^t*expAll*[g'' + 2*g'*u' + g*(u'' + u'^2)]
            double bracket = gpp + 2.0 * gp * u1 + deltaPow * (u2 + u1 * u1);

            d2 += n[i] * tauPow * expAll * bracket;
        }

        return d2;
    }

}
