using LitMath;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Helmholtz;

public static class ResidualHelmholtzGaoBFast
{
    private static readonly double[] b = { 1.244, 0.6826 };
    private static readonly double[] beta = { 0.3696, 0.2962 };
    private static readonly double[] d = { 1, 1 };
    private static readonly double[] epsilon = { 0.4478, 0.44689 };
    private static readonly double[] eta = { -2.8452, -2.8342 };
    private static readonly double[] gamma = { 1.108, 1.313 };
    private static readonly double[] n = { -1.6909858, 0.93739074 };
    private static readonly double[] t = { 4.3315, 4.015 };

    
    private const double b0 = 1.244, b1 = 0.6826;
    private const double beta0 = 0.3696, beta1 = 0.2962;
    private const double eps0 = 0.4478, eps1 = 0.44689;
    private const double eta0 = -2.8452, eta1 = -2.8342;
    private const double gam0 = 1.108, gam1 = 1.313;
    private const double n0 = -1.6909858, n1 = 0.93739074;
    private const double t0 = 4.3315, t1 = 4.015;


    public static double alphaR(double delta, double tau)
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
        // ---------------------------
        // Tau-side shared: log(tau)
        // ---------------------------
        double logTau = Math.Log(tau); // tau must be > 0

        // ---------------------------
        // Term 0: tau side
        // Ftau0 = exp(t0*logTau) * exp(1/denom0)
        // denom0 = b0 + beta0*(gam0-tau)^2
        // ---------------------------
        double dt0 = gam0 - tau;
        double dt0Sq = dt0 * dt0;
        double denom0 = Math.FusedMultiplyAdd(beta0, dt0Sq, b0);
        double expTau0 = Math.Exp(1.0 / denom0);
        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double Ftau0 = tauPow0 * expTau0;

        // ---------------------------
        // Term 1: tau side
        // ---------------------------
        double dt1 = gam1 - tau;
        double dt1Sq = dt1 * dt1;
        double denom1 = Math.FusedMultiplyAdd(beta1, dt1Sq, b1);
        double expTau1 = Math.Exp(1.0 / denom1);
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double Ftau1 = tauPow1 * expTau1;

        // ---------------------------
        // Delta-side
        // expDelta = exp(eta*(delta-eps)^2)
        // multiplier = 1 + 2*eta*delta*(delta-eps)
        // derivative term: n*Ftau*expDelta*multiplier
        // ---------------------------

        // Term 0 delta side
        double dd0 = delta - eps0;
        double dd0Sq = dd0 * dd0;
        double expDelta0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0Sq, 0.0));

        // multiplier0 = 1 + (2*eta0*delta)*dd0
        double k0 = (2.0 * eta0) * delta;
        double mult0 = Math.FusedMultiplyAdd(k0, dd0, 1.0);

        // Term 1 delta side
        double dd1 = delta - eps1;
        double dd1Sq = dd1 * dd1;
        double expDelta1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1Sq, 0.0));

        double k1 = (2.0 * eta1) * delta;
        double mult1 = Math.FusedMultiplyAdd(k1, dd1, 1.0);

        // ---------------------------
        // Combine
        // term = n * Ftau * expDelta * mult
        // Use one FMA for the final sum
        // ---------------------------
        double p0 = (n0 * Ftau0) * expDelta0;   // (n0*Ftau0*expDelta0)
        double p1 = (n1 * Ftau1) * expDelta1;   // (n1*Ftau1*expDelta1)

        // result = p0*mult0 + p1*mult1
        double tail = p1 * mult1;
        return Math.FusedMultiplyAdd(p0, mult0, tail);
    }




    public static double alphaR_dTau(double delta, double tau)
    {
        // tau must be > 0 for Log
        double invTau = 1.0 / tau;
        double logTau = Math.Log(tau);

        // ---------- delta side (depends on i) ----------
        // Fdelta_i = delta * exp(eta_i*(delta-eps_i)^2)
        double dd0 = delta - eps0;
        double dd1 = delta - eps1;

        double expDelta0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0 * dd0, 0.0));
        double expDelta1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1 * dd1, 0.0));

        // Factor delta out later to save 2 multiplies:
        // sum = delta * [ ... ]

        // ---------- term 0 tau side ----------
        double dt0 = gam0 - tau;
        double denom0 = Math.FusedMultiplyAdd(beta0, dt0 * dt0, b0);
        double denom0Sq = denom0 * denom0;

        double expInvDen0 = Math.Exp(1.0 / denom0);
        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double Ftau0 = tauPow0 * expInvDen0;

        // factor0 = t0/tau + 2*beta0*dt0/denom0^2
        double s0 = (2.0 * beta0 * dt0) / denom0Sq;
        double factor0 = Math.FusedMultiplyAdd(t0, invTau, s0);

        // ---------- term 1 tau side ----------
        double dt1 = gam1 - tau;
        double denom1 = Math.FusedMultiplyAdd(beta1, dt1 * dt1, b1);
        double denom1Sq = denom1 * denom1;

        double expInvDen1 = Math.Exp(1.0 / denom1);
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double Ftau1 = tauPow1 * expInvDen1;

        double s1 = (2.0 * beta1 * dt1) / denom1Sq;
        double factor1 = Math.FusedMultiplyAdd(t1, invTau, s1);

        // Combine:
        // term_i = n_i * (delta*expDelta_i) * Ftau_i * factor_i
        // => delta * [ n_i * expDelta_i * Ftau_i * factor_i ]
        double a0 = (n0 * expDelta0) * (Ftau0 * factor0);
        double a1 = (n1 * expDelta1) * (Ftau1 * factor1);

        return delta * (a0 + a1);
    }


    public static double alphaR_dDelta2(double delta, double tau)
    {
        // tau must be > 0 for Log
        double logTau = Math.Log(tau);

        // ---------- tau side (depends on i) ----------
        // Ftau_i = exp(t_i*logTau) * exp(1/denom_i)
        // denom_i = b_i + beta_i*(gamma_i - tau)^2

        // term 0 tau side
        double dt0 = gam0 - tau;
        double denom0 = Math.FusedMultiplyAdd(beta0, dt0 * dt0, b0);
        double expInvDen0 = Math.Exp(1.0 / denom0);
        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double Ftau0 = tauPow0 * expInvDen0;

        // term 1 tau side
        double dt1 = gam1 - tau;
        double denom1 = Math.FusedMultiplyAdd(beta1, dt1 * dt1, b1);
        double expInvDen1 = Math.Exp(1.0 / denom1);
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double Ftau1 = tauPow1 * expInvDen1;

        // ---------- delta side second derivative ----------
        // For d=1:
        // u1 = 2*eta*(delta-eps)
        // u2 = 2*eta
        // E = exp(eta*(delta-eps)^2)
        // Fdelta'' = E * ( 2*u1 + delta*(u2 + u1^2) )
        //
        // We'll do it per term i=0,1

        // term 0 delta side
        double dd0 = delta - eps0;
        double dd0Sq = dd0 * dd0;
        double E0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0Sq, 0.0));

        double u10 = (2.0 * eta0) * dd0;
        double u20 = 2.0 * eta0;
        double bracket0 = Math.FusedMultiplyAdd(delta, (u20 + (u10 * u10)), 2.0 * u10);
        double FdeltaPP0 = E0 * bracket0;

        // term 1 delta side
        double dd1 = delta - eps1;
        double dd1Sq = dd1 * dd1;
        double E1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1Sq, 0.0));

        double u11 = (2.0 * eta1) * dd1;
        double u21 = 2.0 * eta1;
        double bracket1 = Math.FusedMultiplyAdd(delta, (u21 + (u11 * u11)), 2.0 * u11);
        double FdeltaPP1 = E1 * bracket1;

        // Combine: d2 = n0*Ftau0*FdeltaPP0 + n1*Ftau1*FdeltaPP1
        double tail = (n1 * Ftau1) * FdeltaPP1;
        return Math.FusedMultiplyAdd((n0 * Ftau0), FdeltaPP0, tail);
    }




}
