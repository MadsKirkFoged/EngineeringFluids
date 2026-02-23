using System;
using System.Runtime.CompilerServices;
using LitMath; // Lit.Exp(in x, ref y) per LitMath README [3](https://github.com/matthewkolbe/LitMath)[4](https://www.nuget.org/packages/Kolbe.LitMath/)

namespace EngineeringFluids.Helmholtz;

public static class ResidualHelmholtzGaussianSIMD
{
    private static readonly double[] beta = { 1.708, 1.4865, 2.0915, 2.43, 0.488, 1.1, 0.85, 1.14, 945.64, 993.85 };
    private static readonly int[] d = { 1, 1, 1, 2, 2, 1, 3, 3, 1, 1 };
    private static readonly double[] epsilon = { -0.0726, -0.1274, 0.7527, 0.57, 2.2, -0.243, 2.96, 3.02, 0.9574, 0.9576 };
    private static readonly double[] eta = { 0.42776, 0.6424, 0.8175, 0.7995, 0.91, 0.3574, 1.21, 4.14, 22.56, 22.68 };
    private static readonly double[] gamma = { 1.036, 1.2777, 1.083, 1.2906, 0.928, 0.934, 0.919, 1.852, 1.05897, 1.05277 };
    private static readonly double[] n = { 6.212578, -5.7844357, 2.4817542, -2.3739168, 0.01493697, -3.7749264, 0.0006254348, -1.7359e-05, -0.13462033, 0.07749072839 };
    private static readonly double[] t = { 0.655, 1.3, 3.1, 1.4395, 1.623, 0.643, 1.13, 4.5, 1.0, 4.0 };

    private const int N = 10;

    // ---- Thread-local scratch buffers (heap-backed, allocated once per thread) ----
    // We avoid stackalloc because Lit.Exp(in Span<T>, ref Span<T>) + stackalloc triggers CS8350/CS8352. [1](https://stackoverflow.com/questions/61378372/passing-both-stackalloced-spant-and-by-ref-struct-as-arguments)[2](https://blog.walterlv.com/post/cs8350-ref-arguments-combination-is-disallowed)
    [ThreadStatic] private static double[]? _args;    // input args for Lit.Exp
    [ThreadStatic] private static double[]? _tauPow;  // outputs: tau^t
    [ThreadStatic] private static double[]? _expTerm; // outputs: exp(-...)

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void EnsureBuffers()
    {
        if (_args is null || _args.Length < N)
            _args = new double[N];
        if (_tauPow is null || _tauPow.Length < N)
            _tauPow = new double[N];
        if (_expTerm is null || _expTerm.Length < N)
            _expTerm = new double[N];
    }

    // ---- Fast integer powers for delta^d where d ∈ {1,2,3} ----
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double PowDeltaInt(double delta, int di) => di switch
    {
        0 => 1.0,
        1 => delta,
        2 => delta * delta,
        3 => delta * delta * delta,
        _ => Math.Pow(delta, di)
    };

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double PowDeltaIntMinus1(double delta, int di) => di switch
    {
        0 => 0.0,
        1 => 1.0,
        2 => delta,
        3 => delta * delta,
        _ => Math.Pow(delta, di - 1)
    };

    // ---- Compute tauPow[i] = exp(t[i] * logTau) with ONE Lit.Exp call ----
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ComputeTauPow(double tau, double logTau, Span<double> tauPowOut)
    {
        // args[i] = t[i] * logTau
        var argsArr = _args!;
        for (int i = 0; i < N; i++)
        {
            argsArr[i] = Math.FusedMultiplyAdd(t[i], logTau, 0.0);
        }

        // LitMath expects spans; heap-backed spans are fine (no escape issues). [3](https://github.com/matthewkolbe/LitMath)[4](https://www.nuget.org/packages/Kolbe.LitMath/)
        Span<double> inSpan = argsArr.AsSpan(0, N);
        Span<double> outSpan = tauPowOut; // local ref, safe
        Lit.Exp(in inSpan, ref outSpan);  // [3](https://github.com/matthewkolbe/LitMath)[4](https://www.nuget.org/packages/Kolbe.LitMath/)
    }

    // ---- Compute expTerm[i] = exp(-eta*(dd^2) - beta*(tt^2)) with ONE Lit.Exp call ----
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ComputeExpGaussian(double delta, double tau, Span<double> expOut)
    {
        var argsArr = _args!;
        for (int i = 0; i < N; i++)
        {
            double dd = delta - epsilon[i];
            double tt = tau - gamma[i];
            double dd2 = dd * dd;
            double tt2 = tt * tt;

            // arg = -(eta*dd2 + beta*tt2)
            double sum = Math.FusedMultiplyAdd(eta[i], dd2, beta[i] * tt2);
            argsArr[i] = -sum;
        }

        Span<double> inSpan = argsArr.AsSpan(0, N);
        Span<double> outSpan = expOut; // local ref, safe
        Lit.Exp(in inSpan, ref outSpan); // [3](https://github.com/matthewkolbe/LitMath)[4](https://www.nuget.org/packages/Kolbe.LitMath/)
    }

    // ==========================
    // alphaR
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR(double delta, double tau)
    {
        // tau must be > 0 for log (powless)
        double logTau = Math.Log(tau);

        EnsureBuffers();
        Span<double> tauPow = _tauPow!.AsSpan(0, N);
        Span<double> expTerm = _expTerm!.AsSpan(0, N);

        ComputeTauPow(tau, logTau, tauPow);
        ComputeExpGaussian(delta, tau, expTerm);

        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            double deltaPow = PowDeltaInt(delta, d[i]);
            sum += (n[i] * deltaPow) * (tauPow[i] * expTerm[i]);
        }
        return sum;
    }

    // ==========================
    // alphaR_dDelta
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dDelta(double delta, double tau)
    {
        double logTau = Math.Log(tau);

        EnsureBuffers();
        Span<double> tauPow = _tauPow!.AsSpan(0, N);
        Span<double> expTerm = _expTerm!.AsSpan(0, N);

        ComputeTauPow(tau, logTau, tauPow);
        ComputeExpGaussian(delta, tau, expTerm);

        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            int di = d[i];
            double dd = delta - epsilon[i];

            // d/dδ [δ^d * exp(-eta*(δ-eps)^2)] = expTerm * δ^(d-1) * ( d + δ*(-2*eta*dd) )
            double deltaPowDm1 = PowDeltaIntMinus1(delta, di);
            double inner = Math.FusedMultiplyAdd(delta, (-2.0 * eta[i] * dd), di);

            sum += (n[i] * tauPow[i]) * (expTerm[i] * (deltaPowDm1 * inner));
        }
        return sum;
    }

    // ==========================
    // alphaR_dTau
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dTau(double delta, double tau)
    {
        double invTau = 1.0 / tau;
        double logTau = Math.Log(tau);

        EnsureBuffers();
        Span<double> tauPow = _tauPow!.AsSpan(0, N);
        Span<double> expTerm = _expTerm!.AsSpan(0, N);

        ComputeTauPow(tau, logTau, tauPow);
        ComputeExpGaussian(delta, tau, expTerm);

        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            double deltaPow = PowDeltaInt(delta, d[i]);
            double baseTerm = (n[i] * deltaPow) * (tauPow[i] * expTerm[i]);

            double tt = tau - gamma[i];
            double factor = Math.FusedMultiplyAdd(t[i], invTau, (-2.0 * beta[i] * tt));

            sum += baseTerm * factor;
        }
        return sum;
    }

    // ==========================
    // alphaR2_dTau (alternative expression)
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR2_dTau(double delta, double tau)
    {
        double invTau = 1.0 / tau;
        double logTau = Math.Log(tau);

        EnsureBuffers();
        Span<double> tauPow = _tauPow!.AsSpan(0, N);
        Span<double> expTerm = _expTerm!.AsSpan(0, N);

        ComputeTauPow(tau, logTau, tauPow);
        ComputeExpGaussian(delta, tau, expTerm);

        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            double deltaPow = PowDeltaInt(delta, d[i]);

            // NewTerm = t * tau^(t-1) - 2*beta*tau^t*(tau-gamma)
            // tau^(t-1) = tauPow / tau = tauPow * invTau
            double tt = tau - gamma[i];
            double newTerm = Math.FusedMultiplyAdd(t[i], (tauPow[i] * invTau), (-2.0 * beta[i] * tauPow[i] * tt));

            sum += (n[i] * deltaPow) * (newTerm * expTerm[i]);
        }
        return sum;
    }

    // ==========================
    // alphaR_dDelta2
    // ==========================
    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double alphaR_dDelta2(double delta, double tau)
    {
        double logTau = Math.Log(tau);

        EnsureBuffers();
        Span<double> tauPow = _tauPow!.AsSpan(0, N);
        Span<double> expTerm = _expTerm!.AsSpan(0, N);

        ComputeTauPow(tau, logTau, tauPow);
        ComputeExpGaussian(delta, tau, expTerm);

        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            int di = d[i];
            double etai = eta[i];
            double dd = delta - epsilon[i];

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

            // u = -eta*(delta-eps)^2, u' = -2*eta*dd, u'' = -2*eta
            double u1 = (-2.0 * etai) * dd;
            double u2 = -2.0 * etai;

            double u1Sq = u1 * u1;
            double tmp = u2 + u1Sq;

            double bracket = gpp + (2.0 * gp * u1) + (g * tmp);

            sum += (n[i] * tauPow[i]) * (expTerm[i] * bracket);
        }

        return sum;
    }
}