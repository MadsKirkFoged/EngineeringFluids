using System;
using System.Runtime.CompilerServices;

namespace EngineeringFluids.Helmholtz;

public static class ResidualHelmholtzPower
{
    // Coefficients (unrolled as constants to remove loop/array overhead)
    private const double n0 = 0.006132232;
    private const double n1 = 1.7395866;
    private const double n2 = -2.2261792;
    private const double n3 = -0.30127553;
    private const double n4 = 0.08967023;
    private const double n5 = -0.076387037;
    private const double n6 = -0.84063963;
    private const double n7 = -0.27026327;

    // d exponents (integers)
    private const int d0 = 4;
    private const int d1 = 1;
    private const int d2 = 1;
    private const int d3 = 2;
    private const int d4 = 3;
    private const int d5 = 3;
    private const int d6 = 2;
    private const int d7 = 3;

    // t_R exponents (doubles)
    private const double t0 = 1.0;
    private const double t1 = 0.382;
    private const double t2 = 1.0;
    private const double t3 = 1.0;
    private const double t4 = 0.677;
    private const double t5 = 2.915;
    private const double t6 = 3.51;
    private const double t7 = 1.063;

    // l exponents (integers): only 0, 1, 2
    private const int l5 = 2;
    private const int l6 = 2;
    private const int l7 = 1;


    // Every value and derivative in one pass: they all share the same delta powers,
    // exp(-delta^l) dampings and tau^t powers, so these are computed once per state.
    public static ResidualDerivatives Derivatives(double delta, double tau) => Derivatives(delta, tau, Math.Log(tau));

    // logTau = Math.Log(tau), shared by every Helmholtz term of one state
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static ResidualDerivatives Derivatives(double delta, double tau, double logTau)
    {
        // powless tau^t: tau must be > 0
        double invTau = 1.0 / tau;

        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;
        double d5p = d4p * delta; // needed for di+li-1 max 5

        // Exponential damping for l != 0
        double expNegDelta = Math.Exp(-delta); // l = 1
        double expNegDelta2 = Math.Exp(-d2p);   // l = 2

        // tauPow[i] = exp(t[i]*logTau); t0, t2 and t3 are exactly 1, so those are tau itself (as in TauCache)
        double tauPow0 = tau;
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double tauPow2 = tau;
        double tauPow3 = tau;
        double tauPow4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
        double tauPow5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
        double tauPow6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
        double tauPow7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));

        // ---- value: term_i = n * δ^d * τ^t * exp(-δ^l) ----
        double term0 = (n0 * d4p) * tauPow0;
        double term1 = (n1 * d1p) * tauPow1;
        double term2 = (n2 * d1p) * tauPow2;
        double term3 = (n3 * d2p) * tauPow3;
        double term4 = (n4 * d3p) * tauPow4;
        double term5 = (n5 * d3p) * (tauPow5 * expNegDelta2);
        double term6 = (n6 * d2p) * (tauPow6 * expNegDelta2);
        double term7 = (n7 * d3p) * (tauPow7 * expNegDelta);

        double value = term0 + term1 + term2 + term3 + term4 + term5 + term6 + term7;

        // ---- d/dδ ----
        // For l=0: d/dδ [n * δ^d * τ^t] = n * τ^t * d * δ^(d-1)
        double s0 = (n0 * tauPow0) * (4.0 * d3p);  // d=4
        double s1 = (n1 * tauPow1) * 1.0;         // d=1 => δ^0
        double s2 = (n2 * tauPow2) * 1.0;         // d=1
        double s3 = (n3 * tauPow3) * (2.0 * d1p); // d=2 => 2δ
        double s4 = (n4 * tauPow4) * (3.0 * d2p); // d=3 => 3δ^2

        // For l=2: d/dδ [n*τ^t*δ^d*exp(-δ^2)] = n*τ^t*exp(-δ^2) * ( d*δ^(d-1) - 2*δ^(d+1) )
        // term5: d=3 => (3*δ^2 - 2*δ^4)
        double inner5 = Math.FusedMultiplyAdd(-2.0, d4p, 3.0 * d2p);
        double s5 = (n5 * tauPow5) * (expNegDelta2 * inner5);

        // term6: d=2 => (2*δ - 2*δ^3)
        double inner6 = Math.FusedMultiplyAdd(-2.0, d3p, 2.0 * d1p);
        double s6 = (n6 * tauPow6) * (expNegDelta2 * inner6);

        // For l=1: d/dδ [n*τ^t*δ^d*exp(-δ)] = n*τ^t*exp(-δ) * ( d*δ^(d-1) - δ^d )
        // term7: d=3 => (3*δ^2 - δ^3)
        double inner7 = Math.FusedMultiplyAdd(-1.0, d3p, 3.0 * d2p);
        double s7 = (n7 * tauPow7) * (expNegDelta * inner7);

        double dDelta = s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;

        // ---- d/dτ: term * t / τ ----
        double dTau =
            term0 * (t0 * invTau) +
            term1 * (t1 * invTau) +
            term2 * (t2 * invTau) +
            term3 * (t3 * invTau) +
            term4 * (t4 * invTau) +
            term5 * (t5 * invTau) +
            term6 * (t6 * invTau) +
            term7 * (t7 * invTau);

        // ---- d²/dδ² ----
        // l=0 terms: second derivative n*tauPow*d*(d-1)*delta^(d-2)
        double q0 = (n0 * tauPow0) * (12.0 * d2p); // d=4 => 4*3*δ^2
        double q1 = 0.0;                           // d=1 => 0
        double q2 = 0.0;                           // d=1 => 0
        double q3 = (n3 * tauPow3) * 2.0;          // d=2 => 2
        double q4 = (n4 * tauPow4) * (6.0 * d1p);  // d=3 => 6δ

        // l=2 term formula: A*E2*(gpp -4δ gp + g*(4δ^2 -2))
        // term5: d=3 => bracket = 6δ -14δ^3 +4δ^5
        double bracket5 = Math.FusedMultiplyAdd(4.0, d5p, Math.FusedMultiplyAdd(-14.0, d3p, 6.0 * d1p));
        double q5 = (n5 * tauPow5) * (expNegDelta2 * bracket5);

        // term6: d=2 => bracket = 2 -10δ^2 +4δ^4
        double bracket6 = Math.FusedMultiplyAdd(4.0, d4p, (2.0 - 10.0 * d2p));
        double q6 = (n6 * tauPow6) * (expNegDelta2 * bracket6);

        // l=1 term formula: A*E1*(gpp -2gp + g)
        // term7: d=3 => bracket = 6δ -6δ^2 +δ^3
        double bracket7 = Math.FusedMultiplyAdd(1.0, d3p, (6.0 * d1p - 6.0 * d2p));
        double q7 = (n7 * tauPow7) * (expNegDelta * bracket7);

        double dDelta2 = q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7;

        // ---- d²/dτ²: term * t*(t-1) / τ² ----
        double invTau2 = invTau * invTau;
        double dTau2 =
            term0 * (t0 * (t0 - 1.0) * invTau2) +
            term1 * (t1 * (t1 - 1.0) * invTau2) +
            term2 * (t2 * (t2 - 1.0) * invTau2) +
            term3 * (t3 * (t3 - 1.0) * invTau2) +
            term4 * (t4 * (t4 - 1.0) * invTau2) +
            term5 * (t5 * (t5 - 1.0) * invTau2) +
            term6 * (t6 * (t6 - 1.0) * invTau2) +
            term7 * (t7 * (t7 - 1.0) * invTau2);

        // ---- d²/dδdτ: (d/dδ term) * t / τ ----
        double dDeltadTau =
            s0 * (t0 * invTau) +
            s1 * (t1 * invTau) +
            s2 * (t2 * invTau) +
            s3 * (t3 * invTau) +
            s4 * (t4 * invTau) +
            s5 * (t5 * invTau) +
            s6 * (t6 * invTau) +
            s7 * (t7 * invTau);

        return new ResidualDerivatives(value, dDelta, dTau, dDelta2, dTau2, dDeltadTau);
    }

    // Fused first+second delta-derivative: the Newton solver needs both every
    // iteration, and alphaR_dDelta + alphaR_dDelta2 each redundantly evaluate the
    // same 10 Math.Exp calls (2 exp-damping + 8 tauPow) independently. This computes
    // them once and reuses them for both derivatives. [benchmark-guided]
    // Every tau^t power here depends only on tau (i.e. only on temperature), which is
    // fixed for the whole density Newton solve at a given T — only delta changes between
    // iterations. Precomputing these once per solve (instead of once per iteration) removes
    // 1 Log + 8 Exp calls from every Newton step. [benchmark-guided]
    public readonly struct TauCache
    {
        public readonly double p0, p1, p2, p3, p4, p5, p6, p7;

        public TauCache(double tau) : this(tau, Math.Log(tau)) { }

        // logTau is the same value for every residual class's TauCache (they all share
        // the same tau); computing it once in the caller instead of once per class saves
        // 2 redundant Math.Log calls on every single UpdatePT call. [benchmark-guided]
        public TauCache(double tau, double logTau)
        {
            // t0, t2 and t3 are exactly 1.0, so tau^t is just tau - no Exp/Log round-trip
            // needed (also more accurate than exp(log(tau)), not just faster).
            p0 = tau;
            p1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
            p2 = tau;
            p3 = tau;
            p4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
            p5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
            p6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
            p7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));
        }
    }

    // Fused first+second delta-derivative for the density Newton solve, taking a precomputed
    // TauCache instead of recomputing tau^t on every call (see TauCache remarks).
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static void alphaR_dDelta_dDelta2(double delta, in TauCache c, out double dDelta, out double dDelta2)
    {
        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;
        double d5p = d4p * delta;

        double expNegDelta = Math.Exp(-delta);
        double expNegDelta2 = Math.Exp(-d2p);

        // ---- first derivative ----
        double g0 = (n0 * c.p0) * (4.0 * d3p);
        double g1 = (n1 * c.p1) * 1.0;
        double g2 = (n2 * c.p2) * 1.0;
        double g3 = (n3 * c.p3) * (2.0 * d1p);
        double g4 = (n4 * c.p4) * (3.0 * d2p);

        double inner5 = Math.FusedMultiplyAdd(-2.0, d4p, 3.0 * d2p);
        double g5 = (n5 * c.p5) * (expNegDelta2 * inner5);

        double inner6 = Math.FusedMultiplyAdd(-2.0, d3p, 2.0 * d1p);
        double g6 = (n6 * c.p6) * (expNegDelta2 * inner6);

        double inner7 = Math.FusedMultiplyAdd(-1.0, d3p, 3.0 * d2p);
        double g7 = (n7 * c.p7) * (expNegDelta * inner7);

        dDelta = g0 + g1 + g2 + g3 + g4 + g5 + g6 + g7;

        // ---- second derivative (reuses c.p*/expNegDelta*/d*p above) ----
        double s0 = (n0 * c.p0) * (12.0 * d2p);
        double s3 = (n3 * c.p3) * 2.0;
        double s4 = (n4 * c.p4) * (6.0 * d1p);

        double bracket5 = Math.FusedMultiplyAdd(4.0, d5p, Math.FusedMultiplyAdd(-14.0, d3p, 6.0 * d1p));
        double s5 = (n5 * c.p5) * (expNegDelta2 * bracket5);

        double bracket6 = Math.FusedMultiplyAdd(4.0, d4p, (2.0 - 10.0 * d2p));
        double s6 = (n6 * c.p6) * (expNegDelta2 * bracket6);

        double bracket7 = Math.FusedMultiplyAdd(1.0, d3p, (6.0 * d1p - 6.0 * d2p));
        double s7 = (n7 * c.p7) * (expNegDelta * bracket7);

        dDelta2 = s0 + s3 + s4 + s5 + s6 + s7;
    }

    // As above plus alphaR itself - everything a fixed-temperature solve needs (p, dp/drho, ln(phi), d ln(phi)/drho)
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static void alphaR_dDelta_dDelta2(double delta, in TauCache c, out double value, out double dDelta, out double dDelta2)
    {
        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;
        double d5p = d4p * delta;

        double expNegDelta = Math.Exp(-delta);
        double expNegDelta2 = Math.Exp(-d2p);

        value = (n0 * d4p) * c.p0 + (n1 * d1p) * c.p1 + (n2 * d1p) * c.p2 + (n3 * d2p) * c.p3 + (n4 * d3p) * c.p4
              + (n5 * d3p) * (c.p5 * expNegDelta2) + (n6 * d2p) * (c.p6 * expNegDelta2) + (n7 * d3p) * (c.p7 * expNegDelta);

        dDelta = (n0 * c.p0) * (4.0 * d3p) + (n1 * c.p1) + (n2 * c.p2) + (n3 * c.p3) * (2.0 * d1p) + (n4 * c.p4) * (3.0 * d2p)
               + (n5 * c.p5) * (expNegDelta2 * Math.FusedMultiplyAdd(-2.0, d4p, 3.0 * d2p))
               + (n6 * c.p6) * (expNegDelta2 * Math.FusedMultiplyAdd(-2.0, d3p, 2.0 * d1p))
               + (n7 * c.p7) * (expNegDelta * Math.FusedMultiplyAdd(-1.0, d3p, 3.0 * d2p));

        dDelta2 = (n0 * c.p0) * (12.0 * d2p) + (n3 * c.p3) * 2.0 + (n4 * c.p4) * (6.0 * d1p)
                + (n5 * c.p5) * (expNegDelta2 * Math.FusedMultiplyAdd(4.0, d5p, Math.FusedMultiplyAdd(-14.0, d3p, 6.0 * d1p)))
                + (n6 * c.p6) * (expNegDelta2 * Math.FusedMultiplyAdd(4.0, d4p, (2.0 - 10.0 * d2p)))
                + (n7 * c.p7) * (expNegDelta * Math.FusedMultiplyAdd(1.0, d3p, (6.0 * d1p - 6.0 * d2p)));
    }

}
