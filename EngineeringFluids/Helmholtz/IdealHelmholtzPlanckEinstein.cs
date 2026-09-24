using System;
using System.Runtime.CompilerServices;

namespace EngineeringFluids.Helmholtz;

public static class IdealHelmholtzPlanckEinstein
{
    // Original coefficients as scalars (avoid arrays + indexing)
    private const double n0 = 2.224;
    private const double n1 = 3.148;
    private const double n2 = 0.9579;

    // t_PE values (negative)
    private const double t0 = -4.0585856593352405;
    private const double t1 = -9.776605187888352;
    private const double t2 = -17.829667620080876;

    // constants
    private const double Ln2 = 0.693147180559945309417232121458176568; // ln(2)

    // ----------------------------
    // Numerically stable helpers
    // ----------------------------

    /// <summary>
    /// Stable expm1(x) = exp(x) - 1 for small |x|.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Expm1(double x)
    {
        // For very small x, exp(x)-1 loses precision; use a short series.
        // Threshold chosen to keep error tiny while keeping it fast.
        double ax = Math.Abs(x);
        if (ax < 1e-5)
        {
            // 5th-order Taylor: x + x^2/2 + x^3/6 + x^4/24 + x^5/120
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x4 * x;
            return x + 0.5 * x2 + (1.0 / 6.0) * x3 + (1.0 / 24.0) * x4 + (1.0 / 120.0) * x5;
        }
        return Math.Exp(x) - 1.0;
    }

    /// <summary>
    /// Stable log1p(x) = log(1 + x) for small |x|.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Log1p(double x)
    {
        // For |x| small, log(1+x) loses bits due to 1+x rounding; use series.
        double ax = Math.Abs(x);
        if (ax < 1e-4)
        {
            // 5th-order alternating series: x - x^2/2 + x^3/3 - x^4/4 + x^5/5
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x4 * x;
            return x - 0.5 * x2 + (1.0 / 3.0) * x3 - 0.25 * x4 + 0.2 * x5;
        }
        return Math.Log(1.0 + x);
    }

    /// <summary>
    /// Stable log(1 - exp(x)) for x <= 0.
    /// Standard split:
    ///   if x < -ln2  => log1p(-exp(x))
    ///   else         => log(-expm1(x))
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Log1mExp(double x)
    {
        // x should be <= 0 in our usage (t_i negative, tau >= 0)
        if (x < -Ln2)
        {
            // exp(x) is <= 0.5, so -exp(x) is in [-0.5, 0] and log1p is well-behaved.
            return Log1p(-Math.Exp(x));
        }
        else
        {
            // expm1(x) is negative and close to 0; -expm1(x) is small positive.
            return Math.Log(-Expm1(x));
        }
    }

    // ----------------------------
    // Public API
    // ----------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double Alpha0(double delta, double tau)
    {
        // delta is unused (kept for signature compatibility)

        // x_i = t_i * tau (t_i < 0)
        double x0 = t0 * tau;
        double x1 = t1 * tau;
        double x2 = t2 * tau;

        // alpha0 = Σ n_i * log(1 - exp(x_i))   because c=1, d=-1
        return (n0 * Log1mExp(x0)) + (n1 * Log1mExp(x1)) + (n2 * Log1mExp(x2));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static double Alpha0_dTau(double delta, double tau)
    {
        // delta is unused (kept for signature compatibility)

        // Original:
        // Σ n * d * t * exp(x) / (1 + d*exp(x)), with d=-1
        // => Σ n * (-t) * exp(x) / (1 - exp(x))
        //
        // Using expm1:
        // exp(x) = expm1(x) + 1
        // 1 - exp(x) = -(expm1(x))
        //
        // So:
        // (-t)*exp(x)/(1-exp(x)) = (-t)*(expm1+1)/(-expm1) = t*(expm1+1)/expm1
        // => t * exp(x) / expm1(x)
        //
        // This is numerically stable when x ~ 0- because expm1(x) keeps precision.

        double x0 = t0 * tau;
        double x1 = t1 * tau;
        double x2 = t2 * tau;

        double em10 = Expm1(x0);
        double em11 = Expm1(x1);
        double em12 = Expm1(x2);

        double e0 = em10 + 1.0;
        double e1 = em11 + 1.0;
        double e2 = em12 + 1.0;

        double s0 = n0 * t0 * (e0 / em10);
        double s1 = n1 * t1 * (e1 / em11);
        double s2 = n2 * t2 * (e2 / em12);

        return s0 + s1 + s2;
    }
}