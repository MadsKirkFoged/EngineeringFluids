using System;
using System.Runtime.CompilerServices;

namespace EngineeringFluids.Helmholtz;

public static class ResidualHelmholtzGaussian
{
    // Keep original coefficients (same numbers), but also keep an int[] version for d
    // so we can avoid (int)d[i] casts in hot loops.
    private static readonly double[] beta =
    {
        1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85
    };

    private static readonly int[] dInt =
    {
        1, 1, 1, 2, 2, 1, 3, 3, 1, 1
    };

    private static readonly double[] epsilon =
    {
        -0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576
    };

    private static readonly double[] eta =
    {
        0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68
    };

    private static readonly double[] gamma =
    {
        1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277
    };

    private static readonly double[] n =
    {
        6.212578, -5.7844357, 2.4817542, -2.3739168, 0.01493697, -3.7749264,
        0.0006254348, -1.7359e-05, -0.13462033, 0.07749072839
    };

    private static readonly double[] t =
    {
        0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13, 4.5, 1.0, 4.0
    };

    private const int N = 10;

    // --------------------------
    // Small helper: delta^d where d is 1,2,3 (fast path). Fallback to Math.Pow if ever needed.
    // --------------------------
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double PowIntDelta(double delta, int di)
    {
        return di switch
        {
            0 => 1.0,
            1 => delta,
            2 => delta * delta,
            3 => delta * delta * delta,
            _ => Math.Pow(delta, di) // should never hit with current dInt values
        };
    }

    // --------------------------
    // Small helper: delta^(d-1) for d=1,2,3 (fast path)
    // --------------------------
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double PowIntDeltaMinus1(double delta, int di)
    {
        return di switch
        {
            0 => 0.0,             // derivative factor will zero out anyway
            1 => 1.0,
            2 => delta,
            3 => delta * delta,
            _ => di * Math.Pow(delta, di - 1) / di // fallback (avoids duplication)
        };
    }

    // --------------------------
    // Helper: tau^t via powless exp(t*logTau). tau > 0 required.
    // --------------------------
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double TauPow(double ti, double logTau)
        => Math.Exp(Math.FusedMultiplyAdd(ti, logTau, 0.0));

    // --------------------------
    // Helper: exp( -eta*(delta-eps)^2 - beta*(tau-gam)^2 )
    // Uses FMA for (eta*dd2 + beta*tt2) then negates.
    // --------------------------
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double ExpGaussian(double etai, double dd, double betai, double tt)
    {
        double dd2 = dd * dd;
        double tt2 = tt * tt;
        double sum = Math.FusedMultiplyAdd(etai, dd2, betai * tt2); // etai*dd2 + betai*tt2
        return Math.Exp(-sum);
    }

    // ==========================
    // Derivatives
    // Every value and derivative in one pass over the N terms: they all share the same
    // tau^t and exp(-eta*(delta-eps)^2 - beta*(tau-gamma)^2) factors, so each term's
    // Exp/Log work is done once per state. Value/dDelta/dTau/dDelta2 keep the exact
    // arithmetic of the former separate alphaR/alphaR_dDelta/alphaR_dTau/alphaR_dDelta2.
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static ResidualDerivatives Derivatives(double delta, double tau)
    {
        // tau must be > 0 for log
        double invTau = 1.0 / tau;
        double logTau = Math.Log(tau);

        // Get refs for Unsafe.Add access (bounds-check elimination)
        ref double betaRef = ref beta[0];
        ref int dRef = ref dInt[0];
        ref double epsRef = ref epsilon[0];
        ref double etaRef = ref eta[0];
        ref double gamRef = ref gamma[0];
        ref double nRef = ref n[0];
        ref double tRef = ref t[0];

        double sumValue = 0.0;
        double sumDelta = 0.0;
        double sumTau = 0.0;
        double sumDelta2 = 0.0;
        double sumTau2 = 0.0;
        double sumDeltaTau = 0.0;

        for (int i = 0; i < N; i++)
        {
            double betai = Unsafe.Add(ref betaRef, i);
            int di = Unsafe.Add(ref dRef, i);
            double epsi = Unsafe.Add(ref epsRef, i);
            double etai = Unsafe.Add(ref etaRef, i);
            double gami = Unsafe.Add(ref gamRef, i);
            double ni = Unsafe.Add(ref nRef, i);
            double ti = Unsafe.Add(ref tRef, i);

            double deltaPow = PowIntDelta(delta, di);
            double tauPow = TauPow(ti, logTau);

            double dd = delta - epsi;
            double tt = tau - gami;

            // expTerm = exp( -eta*dd^2 - beta*tt^2 )
            double expTerm = ExpGaussian(etai, dd, betai, tt);

            // ---- value: term = n * delta^d * tau^t * expTerm ----
            double term = (ni * deltaPow) * (tauPow * expTerm);
            sumValue += term;

            // ---- d/dδ ----
            // d/dδ [δ^d * exp(-eta*(δ-eps)^2)] = expTerm * δ^(d-1) * ( d + δ*(-2*eta*(δ-eps)) )
            double deltaPowDm1 = PowIntDeltaMinus1(delta, di);

            // inner = d + delta * (-2*eta*dd)
            double inner = Math.FusedMultiplyAdd(delta, (-2.0 * etai * dd), di);

            double termDelta = (ni * tauPow) * (expTerm * (deltaPowDm1 * inner));
            sumDelta += termDelta;

            // ---- d/dτ: term * ( t/tau - 2*beta*(tau-gamma) ) ----
            double factor = Math.FusedMultiplyAdd(ti, invTau, (-2.0 * betai * tt));
            sumTau += term * factor;

            // ---- d²/dδ² ----
            // g = delta^d, gp = d*delta^(d-1), gpp = d*(d-1)*delta^(d-2)
            // d in {1,2,3} so do it branchlessly-ish via switch
            double g, gp, gpp;
            switch (di)
            {
                case 1:
                    g = delta;
                    gp = 1.0;
                    gpp = 0.0;
                    break;

                case 2:
                    g = delta * delta;
                    gp = 2.0 * delta;
                    gpp = 2.0;
                    break;

                case 3:
                    double d2 = delta * delta;
                    g = d2 * delta;
                    gp = 3.0 * d2;
                    gpp = 6.0 * delta;
                    break;

                default:
                    // Should never happen with current dInt, but keep correctness.
                    g = Math.Pow(delta, di);
                    gp = (di == 0) ? 0.0 : di * Math.Pow(delta, di - 1);
                    gpp = (di <= 1) ? 0.0 : di * (di - 1) * Math.Pow(delta, di - 2);
                    break;
            }

            // u(delta) = -eta*(delta-eps)^2, u' = -2*eta*(delta-eps), u'' = -2*eta
            double u1 = (-2.0 * etai) * dd;
            double u2 = -2.0 * etai;

            // bracket = g'' + 2*g'*u' + g*(u'' + u'^2)
            double u1Sq = u1 * u1;
            double tmp = u2 + u1Sq;
            double bracket = gpp + (2.0 * gp * u1) + (g * tmp);

            sumDelta2 += (ni * tauPow) * (expTerm * bracket);

            // ---- d²/dτ²: term * ( factor^2 + d(factor)/dτ ), d(factor)/dτ = -t/tau^2 - 2*beta ----
            sumTau2 += term * (factor * factor - ti * invTau * invTau - 2.0 * betai);

            // ---- d²/dδdτ: (d/dδ term) * factor ----
            sumDeltaTau += termDelta * factor;
        }

        return new ResidualDerivatives(sumValue, sumDelta, sumTau, sumDelta2, sumTau2, sumDeltaTau);
    }

    [System.Runtime.CompilerServices.InlineArray(N)]
    public struct Buffer10
    {
        private double _e0;
    }

    // Every term's tau^t and exp(-beta*(tau-gamma)^2) factor depends only on tau (i.e. only
    // on temperature), which is fixed for the whole density Newton solve at a given T - only
    // delta changes between iterations. Precomputing both once per solve (instead of once per
    // iteration) removes 1 Log + 10 Exp calls from every Newton step; the combined
    // exp(-eta*dd^2-beta*tt^2) call becomes exp(-eta*dd^2) * (precomputed exp(-beta*tt^2)).
    //
    // The two precomputed factors are themselves only ever used multiplied together
    // (tau^t * exp(-beta*(tau-gamma)^2)), so exp(a)*exp(b) = exp(a+b) lets them be folded
    // into ONE Exp call per term instead of two, halving this cache's setup cost (10 Exp
    // instead of 20). This setup runs once per UpdatePT call regardless of how many Newton
    // iterations follow, so it benefits every phase (gas/liquid/supercritical) equally -
    // for the cheap gas/liquid cases (2-4 iterations) this setup was already the majority
    // of the call's cost. [benchmark-guided]
    public readonly struct TauCache
    {
        public readonly Buffer10 Ftau;

        public TauCache(double tau) : this(tau, Math.Log(tau)) { }

        public TauCache(double tau, double logTau)
        {
            ref double betaRef = ref beta[0];
            ref double gamRef = ref gamma[0];
            ref double tRef = ref t[0];

            for (int i = 0; i < N; i++)
            {
                double betai = Unsafe.Add(ref betaRef, i);
                double gami = Unsafe.Add(ref gamRef, i);
                double ti = Unsafe.Add(ref tRef, i);

                double tt = tau - gami;
                double exponent = Math.FusedMultiplyAdd(ti, logTau, -(betai * tt * tt));
                Ftau[i] = Math.Exp(exponent);
            }
        }
    }

    // ==========================
    // alphaR_dDelta_dDelta2 (precomputed TauCache overload)
    // Fused first+second delta-derivative for the density Newton solve that only recomputes
    // the delta-dependent exponential each iteration (see TauCache remarks).
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static void alphaR_dDelta_dDelta2(double delta, in TauCache cache, out double dDelta, out double dDelta2)
    {
        ref int dRef = ref dInt[0];
        ref double epsRef = ref epsilon[0];
        ref double etaRef = ref eta[0];
        ref double nRef = ref n[0];

        double sum1 = 0.0;
        double sum2 = 0.0;

        for (int i = 0; i < N; i++)
        {
            int di = Unsafe.Add(ref dRef, i);
            double epsi = Unsafe.Add(ref epsRef, i);
            double etai = Unsafe.Add(ref etaRef, i);
            double ni = Unsafe.Add(ref nRef, i);

            double dd = delta - epsi;

            // cache.Ftau[i] = tau^t * exp(-beta*(tau-gamma)^2), precomputed and constant
            // across this Newton solve; only the delta-dependent factor is fresh here.
            double expEtaDelta = Math.Exp(-(etai * dd * dd));
            double niTauPowExp = ni * cache.Ftau[i] * expEtaDelta;

            // --- first derivative wrt delta ---
            double deltaPowDm1 = PowIntDeltaMinus1(delta, di);
            double inner = Math.FusedMultiplyAdd(delta, (-2.0 * etai * dd), di);
            sum1 += niTauPowExp * (deltaPowDm1 * inner);

            // --- second derivative wrt delta ---
            double g, gp, gpp;
            switch (di)
            {
                case 1:
                    g = delta;
                    gp = 1.0;
                    gpp = 0.0;
                    break;

                case 2:
                    g = delta * delta;
                    gp = 2.0 * delta;
                    gpp = 2.0;
                    break;

                case 3:
                    double d2 = delta * delta;
                    g = d2 * delta;
                    gp = 3.0 * d2;
                    gpp = 6.0 * delta;
                    break;

                default:
                    g = Math.Pow(delta, di);
                    gp = (di == 0) ? 0.0 : di * Math.Pow(delta, di - 1);
                    gpp = (di <= 1) ? 0.0 : di * (di - 1) * Math.Pow(delta, di - 2);
                    break;
            }

            double u1 = (-2.0 * etai) * dd;
            double u2 = -2.0 * etai;
            double u1Sq = u1 * u1;
            double tmp = u2 + u1Sq;
            double bracket = gpp + (2.0 * gp * u1) + (g * tmp);

            sum2 += niTauPowExp * bracket;
        }

        dDelta = sum1;
        dDelta2 = sum2;
    }

}
