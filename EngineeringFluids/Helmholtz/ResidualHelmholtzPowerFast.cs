using System;
using System.Runtime.CompilerServices;

namespace EngineeringFluids.Helmholtz;

public static class ResidualHelmholtzPowerFast
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

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR(double delta, double tau)
    {
        // powless tau^t: tau must be > 0
        double logTau = Math.Log(tau);

        // Precompute delta powers up to 5 (needed in derivatives)
        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;

        // Exponential damping for l != 0
        double expNegDelta = Math.Exp(-delta); // l = 1
        double expNegDelta2 = Math.Exp(-d2p);   // l = 2

        // tauPow[i] = exp(t[i]*logTau)
        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double tauPow2 = Math.Exp(Math.FusedMultiplyAdd(t2, logTau, 0.0));
        double tauPow3 = Math.Exp(Math.FusedMultiplyAdd(t3, logTau, 0.0));
        double tauPow4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
        double tauPow5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
        double tauPow6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
        double tauPow7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));

        // Unrolled sum
        double sum =
            (n0 * d4p) * tauPow0 +
            (n1 * d1p) * tauPow1 +
            (n2 * d1p) * tauPow2 +
            (n3 * d2p) * tauPow3 +
            (n4 * d3p) * tauPow4 +
            (n5 * d3p) * (tauPow5 * expNegDelta2) +
            (n6 * d2p) * (tauPow6 * expNegDelta2) +
            (n7 * d3p) * (tauPow7 * expNegDelta);

        return sum;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dDelta(double delta, double tau)
    {
        double logTau = Math.Log(tau);

        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;
        double d5p = d4p * delta; // needed for di+li-1 max 5

        double expNegDelta = Math.Exp(-delta); // l=1
        double expNegDelta2 = Math.Exp(-d2p);   // l=2

        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double tauPow2 = Math.Exp(Math.FusedMultiplyAdd(t2, logTau, 0.0));
        double tauPow3 = Math.Exp(Math.FusedMultiplyAdd(t3, logTau, 0.0));
        double tauPow4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
        double tauPow5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
        double tauPow6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
        double tauPow7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));

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

        return s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dTau(double delta, double tau)
    {
        double invTau = 1.0 / tau;
        double logTau = Math.Log(tau);

        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;

        double expNegDelta = Math.Exp(-delta);
        double expNegDelta2 = Math.Exp(-d2p);

        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double tauPow2 = Math.Exp(Math.FusedMultiplyAdd(t2, logTau, 0.0));
        double tauPow3 = Math.Exp(Math.FusedMultiplyAdd(t3, logTau, 0.0));
        double tauPow4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
        double tauPow5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
        double tauPow6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
        double tauPow7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));

        // base term_i = n * δ^d * τ^t * E
        double term0 = (n0 * d4p) * tauPow0;
        double term1 = (n1 * d1p) * tauPow1;
        double term2 = (n2 * d1p) * tauPow2;
        double term3 = (n3 * d2p) * tauPow3;
        double term4 = (n4 * d3p) * tauPow4;
        double term5 = (n5 * d3p) * (tauPow5 * expNegDelta2);
        double term6 = (n6 * d2p) * (tauPow6 * expNegDelta2);
        double term7 = (n7 * d3p) * (tauPow7 * expNegDelta);

        // d/dτ term = term * t / τ
        double sum =
            term0 * (t0 * invTau) +
            term1 * (t1 * invTau) +
            term2 * (t2 * invTau) +
            term3 * (t3 * invTau) +
            term4 * (t4 * invTau) +
            term5 * (t5 * invTau) +
            term6 * (t6 * invTau) +
            term7 * (t7 * invTau);

        return sum;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR2_dTau(double delta, double tau)
    {
        // This is mathematically identical to alphaR_dTauFast:
        // n*δ^d*t*τ^(t-1)*E = (n*δ^d*τ^t*E) * t/τ
        return alphaR_dTau(delta, tau);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dDelta2(double delta, double tau)
    {
        double logTau = Math.Log(tau);

        double d1p = delta;
        double d2p = delta * delta;
        double d3p = d2p * delta;
        double d4p = d2p * d2p;
        double d5p = d4p * delta;

        double expNegDelta = Math.Exp(-delta);
        double expNegDelta2 = Math.Exp(-d2p);

        double tauPow0 = Math.Exp(Math.FusedMultiplyAdd(t0, logTau, 0.0));
        double tauPow1 = Math.Exp(Math.FusedMultiplyAdd(t1, logTau, 0.0));
        double tauPow2 = Math.Exp(Math.FusedMultiplyAdd(t2, logTau, 0.0));
        double tauPow3 = Math.Exp(Math.FusedMultiplyAdd(t3, logTau, 0.0));
        double tauPow4 = Math.Exp(Math.FusedMultiplyAdd(t4, logTau, 0.0));
        double tauPow5 = Math.Exp(Math.FusedMultiplyAdd(t5, logTau, 0.0));
        double tauPow6 = Math.Exp(Math.FusedMultiplyAdd(t6, logTau, 0.0));
        double tauPow7 = Math.Exp(Math.FusedMultiplyAdd(t7, logTau, 0.0));

        // l=0 terms: second derivative n*tauPow*d*(d-1)*delta^(d-2)
        double s0 = (n0 * tauPow0) * (12.0 * d2p); // d=4 => 4*3*δ^2
        double s1 = 0.0;                           // d=1 => 0
        double s2 = 0.0;                           // d=1 => 0
        double s3 = (n3 * tauPow3) * 2.0;          // d=2 => 2
        double s4 = (n4 * tauPow4) * (6.0 * d1p);  // d=3 => 6δ

        // l=2 term formula: A*E2*(gpp -4δ gp + g*(4δ^2 -2))
        // term5: d=3 => bracket = 6δ -14δ^3 +4δ^5
        double bracket5 = Math.FusedMultiplyAdd(4.0, d5p, Math.FusedMultiplyAdd(-14.0, d3p, 6.0 * d1p));
        double s5 = (n5 * tauPow5) * (expNegDelta2 * bracket5);

        // term6: d=2 => bracket = 2 -10δ^2 +4δ^4
        double bracket6 = Math.FusedMultiplyAdd(4.0, d4p, (2.0 - 10.0 * d2p));
        double s6 = (n6 * tauPow6) * (expNegDelta2 * bracket6);

        // l=1 term formula: A*E1*(gpp -2gp + g)
        // term7: d=3 => bracket = 6δ -6δ^2 +δ^3
        double bracket7 = Math.FusedMultiplyAdd(1.0, d3p, (6.0 * d1p - 6.0 * d2p));
        double s7 = (n7 * tauPow7) * (expNegDelta * bracket7);

        return s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;
    }
}