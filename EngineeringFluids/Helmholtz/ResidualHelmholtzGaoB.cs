using LitMath;
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




    private const double b0 = 1.244, b1 = 0.6826;
    private const double beta0 = 0.3696, beta1 = 0.2962;
    private const double eps0 = 0.4478, eps1 = 0.44689;
    private const double eta0 = -2.8452, eta1 = -2.8342;
    private const double gam0 = 1.108, gam1 = 1.313;
    private const double n0 = -1.6909858, n1 = 0.93739074;
    private const double t0 = 4.3315, t1 = 4.015;


    public static double alphaRFast(double delta, double tau)
    {
        // ---------------------------
        // Delta side
        // ---------------------------
        double dd0 = delta - eps0;
        double dd1 = delta - eps1;

        double dd0Sq = dd0 * dd0;
        double dd1Sq = dd1 * dd1;

        // expDelta = exp(eta * dd^2)
        // Use FMA form (eta*ddSq + 0) to encourage fused multiply-add instruction selection where applicable.
        double expDelta0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0Sq, 0.0));
        double expDelta1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1Sq, 0.0));

        // We'll factor delta out later: final = delta * ( ... )
        // So keep expDelta0/1 as-is.

        // ---------------------------
        // Tau side (Powless)
        // ---------------------------
        double dt0 = gam0 - tau;
        double dt1 = gam1 - tau;

        double dt0Sq = dt0 * dt0;
        double dt1Sq = dt1 * dt1;

        // denom = b + beta * dt^2  (FMA is ideal here)
        double denom0 = Math.FusedMultiplyAdd(beta0, dt0Sq, b0);
        double denom1 = Math.FusedMultiplyAdd(beta1, dt1Sq, b1);

        // expTau = exp(1/denom)
        double expTau0 = Math.Exp(1.0 / denom0);
        double expTau1 = Math.Exp(1.0 / denom1);

        // tau^t = exp(t * log(tau)) with one log
        double logTau = Math.Log(tau);

        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));

        // Ftau = tau^t * exp(1/denom)
        double ftau0 = tauPow0 * expTau0;
        double ftau1 = tauPow1 * expTau1;

        // ---------------------------
        // Final combine
        // Original: (n0*ftau0*(delta*expDelta0)) + (n1*ftau1*(delta*expDelta1))
        // Factor delta out:
        //   delta * [ (n0*ftau0*expDelta0) + (n1*ftau1*expDelta1) ]
        //
        // Use FMA for the sum:
        //   sum = FMA(n0*ftau0, expDelta0, (n1*ftau1*expDelta1))
        // ---------------------------
        double a0 = (n0 * ftau0);
        double a1 = (n1 * ftau1);

        double tail = a1 * expDelta1;
        double sum = Math.FusedMultiplyAdd(a0, expDelta0, tail);

        return delta * sum;
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
