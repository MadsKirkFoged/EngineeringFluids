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

    
    private const double b0 = 1.244, b1 = 0.6826;
    private const double beta0 = 0.3696, beta1 = 0.2962;
    private const double eps0 = 0.4478, eps1 = 0.44689;
    private const double eta0 = -2.8452, eta1 = -2.8342;
    private const double gam0 = 1.108, gam1 = 1.313;
    private const double n0 = -1.6909858, n1 = 0.93739074;
    private const double t0 = 4.3315, t1 = 4.015;

    // Every value and derivative in one pass: they share the same tau-side (tau^t,
    // exp(1/denom)) and delta-side (exp(eta*dd^2)) exponentials, so these are computed once
    // per state. Value/dDelta/dTau/dDelta2 keep the exact arithmetic of the former separate
    // alphaR/alphaR_dDelta/alphaR_dTau/alphaR_dDelta2 functions.
    //
    // term_i = n_i * delta * tau^t_i * exp(eta_i*(delta-eps_i)^2 + 1/(b_i + beta_i*(gamma_i-tau)^2))
    public static ResidualDerivatives Derivatives(double delta, double tau) => Derivatives(delta, tau, Math.Log(tau));

    // logTau = Math.Log(tau), shared by every Helmholtz term of one state
    public static ResidualDerivatives Derivatives(double delta, double tau, double logTau)
    {
        // tau must be > 0 for Log
        double invTau = 1.0 / tau;

        // ---------------------------
        // Delta side: expDelta = exp(eta * dd^2)
        // ---------------------------
        double dd0 = delta - eps0;
        double dd1 = delta - eps1;

        double dd0Sq = dd0 * dd0;
        double dd1Sq = dd1 * dd1;

        double expDelta0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0Sq, 0.0));
        double expDelta1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1Sq, 0.0));

        // ---------------------------
        // Tau side: Ftau = tau^t * exp(1/denom), denom = b + beta * (gam - tau)^2
        // ---------------------------
        double dt0 = gam0 - tau;
        double dt1 = gam1 - tau;

        double dt0Sq = dt0 * dt0;
        double dt1Sq = dt1 * dt1;

        double denom0 = Math.FusedMultiplyAdd(beta0, dt0Sq, b0);
        double denom1 = Math.FusedMultiplyAdd(beta1, dt1Sq, b1);

        // tau^t * exp(1/denom) in one Exp each (as in TauCache)
        double Ftau0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 1.0 / denom0));
        double Ftau1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 1.0 / denom1));

        // ---- value = delta * [ (n0*Ftau0*expDelta0) + (n1*Ftau1*expDelta1) ] ----
        double a0 = (n0 * Ftau0);
        double a1 = (n1 * Ftau1);
        double value = delta * Math.FusedMultiplyAdd(a0, expDelta0, a1 * expDelta1);

        // ---- d/dδ: n*Ftau*expDelta*(1 + 2*eta*delta*(delta-eps)) ----
        double k0 = (2.0 * eta0) * delta;
        double mult0 = Math.FusedMultiplyAdd(k0, dd0, 1.0);

        double k1 = (2.0 * eta1) * delta;
        double mult1 = Math.FusedMultiplyAdd(k1, dd1, 1.0);

        double p0 = (n0 * Ftau0) * expDelta0;
        double p1 = (n1 * Ftau1) * expDelta1;

        double dDelta = Math.FusedMultiplyAdd(p0, mult0, p1 * mult1);

        // ---- d/dτ: term * factor, factor = t/tau + 2*beta*(gam-tau)/denom^2 ----
        double denom0Sq = denom0 * denom0;
        double denom1Sq = denom1 * denom1;

        double factor0 = Math.FusedMultiplyAdd(t0, invTau, (2.0 * beta0 * dt0) / denom0Sq);
        double factor1 = Math.FusedMultiplyAdd(t1, invTau, (2.0 * beta1 * dt1) / denom1Sq);

        double dTau = delta * (((n0 * expDelta0) * (Ftau0 * factor0)) + ((n1 * expDelta1) * (Ftau1 * factor1)));

        // ---- d²/dδ²: n*Ftau*expDelta*(2*u1 + delta*(u2 + u1^2)), u1 = 2*eta*dd, u2 = 2*eta ----
        double u10 = (2.0 * eta0) * dd0;
        double u20 = 2.0 * eta0;
        double bracket0 = Math.FusedMultiplyAdd(delta, (u20 + (u10 * u10)), 2.0 * u10);
        double FdeltaPP0 = expDelta0 * bracket0;

        double u11 = (2.0 * eta1) * dd1;
        double u21 = 2.0 * eta1;
        double bracket1 = Math.FusedMultiplyAdd(delta, (u21 + (u11 * u11)), 2.0 * u11);
        double FdeltaPP1 = expDelta1 * bracket1;

        double dDelta2 = Math.FusedMultiplyAdd((n0 * Ftau0), FdeltaPP0, (n1 * Ftau1) * FdeltaPP1);

        // ---- d²/dτ²: term * (factor^2 + d(factor)/dτ) ----
        // d(factor)/dτ = -t/tau^2 + (8*beta^2*(gam-tau)^2 - 2*beta*denom) / denom^3
        double invTau2 = invTau * invTau;
        double dFactor0 = -t0 * invTau2 + (8.0 * beta0 * beta0 * dt0Sq - 2.0 * beta0 * denom0) / (denom0Sq * denom0);
        double dFactor1 = -t1 * invTau2 + (8.0 * beta1 * beta1 * dt1Sq - 2.0 * beta1 * denom1) / (denom1Sq * denom1);

        double dTau2 = delta * (((n0 * expDelta0) * (Ftau0 * (factor0 * factor0 + dFactor0)))
                              + ((n1 * expDelta1) * (Ftau1 * (factor1 * factor1 + dFactor1))));

        // ---- d²/dδdτ: (d/dδ term) * factor ----
        double dDeltadTau = (p0 * mult0) * factor0 + (p1 * mult1) * factor1;

        return new ResidualDerivatives(value, dDelta, dTau, dDelta2, dTau2, dDeltadTau);
    }

    // Ftau0/Ftau1 depend only on tau (i.e. only on temperature), which is fixed for the
    // whole density Newton solve at a given T - only delta changes between iterations.
    // Precomputing them once per solve (instead of once per iteration) removes 1 Log + 4 Exp
    // calls from every Newton step. Each Ftau is tau^t * exp(1/denom) - since exp(a)*exp(b)
    // = exp(a+b), that product is folded into a single Exp call per term instead of two,
    // halving this cache's setup cost too. [benchmark-guided]
    public readonly struct TauCache
    {
        public readonly double Ftau0, Ftau1;

        public TauCache(double tau) : this(tau, Math.Log(tau)) { }

        public TauCache(double tau, double logTau)
        {
            double dt0 = gam0 - tau;
            double denom0 = Math.FusedMultiplyAdd(beta0, dt0 * dt0, b0);
            Ftau0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 1.0 / denom0));

            double dt1 = gam1 - tau;
            double denom1 = Math.FusedMultiplyAdd(beta1, dt1 * dt1, b1);
            Ftau1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 1.0 / denom1));
        }
    }

    // Fused first+second delta-derivative for the density Newton solve, taking a precomputed
    // TauCache instead of recomputing the tau-side terms on every call (see TauCache remarks).
    public static void alphaR_dDelta_dDelta2(double delta, in TauCache c, out double dDelta, out double dDelta2)
        => alphaR_dDelta_dDelta2(delta, in c, out _, out dDelta, out dDelta2);

    // As above plus alphaR itself - everything a fixed-temperature solve needs (p, dp/drho, ln(phi), d ln(phi)/drho)
    public static void alphaR_dDelta_dDelta2(double delta, in TauCache c, out double value, out double dDelta, out double dDelta2)
    {
        // term 0 delta side
        double dd0 = delta - eps0;
        double dd0Sq = dd0 * dd0;
        double expDelta0 = Math.Exp(Math.FusedMultiplyAdd(eta0, dd0Sq, 0.0));

        // term 1 delta side
        double dd1 = delta - eps1;
        double dd1Sq = dd1 * dd1;
        double expDelta1 = Math.Exp(Math.FusedMultiplyAdd(eta1, dd1Sq, 0.0));

        // ---- first derivative ----
        double k0 = (2.0 * eta0) * delta;
        double mult0 = Math.FusedMultiplyAdd(k0, dd0, 1.0);

        double k1 = (2.0 * eta1) * delta;
        double mult1 = Math.FusedMultiplyAdd(k1, dd1, 1.0);

        double p0 = (n0 * c.Ftau0) * expDelta0;
        double p1 = (n1 * c.Ftau1) * expDelta1;

        dDelta = Math.FusedMultiplyAdd(p0, mult0, p1 * mult1);
        value = delta * (p0 + p1);

        // ---- second derivative (reuses expDelta0/1, dd0/1) ----
        double u10 = (2.0 * eta0) * dd0;
        double u20 = 2.0 * eta0;
        double bracket0 = Math.FusedMultiplyAdd(delta, (u20 + (u10 * u10)), 2.0 * u10);
        double FdeltaPP0 = expDelta0 * bracket0;

        double u11 = (2.0 * eta1) * dd1;
        double u21 = 2.0 * eta1;
        double bracket1 = Math.FusedMultiplyAdd(delta, (u21 + (u11 * u11)), 2.0 * u11);
        double FdeltaPP1 = expDelta1 * bracket1;

        dDelta2 = Math.FusedMultiplyAdd(c.Ftau0 * n0, FdeltaPP0, (n1 * c.Ftau1) * FdeltaPP1);
    }

}
