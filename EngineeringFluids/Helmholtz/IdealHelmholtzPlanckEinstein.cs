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
    /// Stable expm1(x) = exp(x) - 1 for small |x|; <paramref name="expX"/> is Math.Exp(x), already computed by the caller.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Expm1(double x, double expX)
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
        return expX - 1.0;
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
    /// Stable log(1 - exp(x)) for x &lt;= 0.
    /// Standard split:
    ///   if x &lt; -ln2  =&gt; log1p(-exp(x))
    ///   else         =&gt; log(-expm1(x))
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Log1mExp(double x, double expX)
    {
        // x should be <= 0 in our usage (t_i negative, tau >= 0)
        if (x < -Ln2)
        {
            // exp(x) is <= 0.5, so -exp(x) is in [-0.5, 0] and log1p is well-behaved.
            return Log1p(-expX);
        }
        else
        {
            // expm1(x) is negative and close to 0; -expm1(x) is small positive.
            return Math.Log(-Expm1(x, expX));
        }
    }

    // ----------------------------
    // Public API
    // ----------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static IdealDerivatives Derivatives(double delta, double tau)
    {
        // delta is unused (kept for signature compatibility)

        // x_i = t_i * tau (t_i < 0)
        double x0 = t0 * tau;
        double x1 = t1 * tau;
        double x2 = t2 * tau;

        // exp(x_i) once per term, shared by the value (Log1mExp) and the derivatives (Expm1) - bit-identical to
        // calling both separately, which evaluated the same Math.Exp twice
        double ex0 = Math.Exp(x0), ex1 = Math.Exp(x1), ex2 = Math.Exp(x2);

        // alpha0 = Σ n_i * log(1 - exp(x_i))   because c=1, d=-1
        double value = (n0 * Log1mExp(x0, ex0)) + (n1 * Log1mExp(x1, ex1)) + (n2 * Log1mExp(x2, ex2));

        // d/dtau:
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
        double em10 = Expm1(x0, ex0);
        double em11 = Expm1(x1, ex1);
        double em12 = Expm1(x2, ex2);

        double e0 = em10 + 1.0;
        double e1 = em11 + 1.0;
        double e2 = em12 + 1.0;

        double s0 = n0 * t0 * (e0 / em10);
        double s1 = n1 * t1 * (e1 / em11);
        double s2 = n2 * t2 * (e2 / em12);

        // d2/dtau2 of n*t*exp(x)/expm1(x) = -n * t^2 * exp(x) / expm1(x)^2
        double q0 = -(n0 * t0 * t0) * (e0 / (em10 * em10));
        double q1 = -(n1 * t1 * t1) * (e1 / (em11 * em11));
        double q2 = -(n2 * t2 * t2) * (e2 / (em12 * em12));

        return new IdealDerivatives(value, s0 + s1 + s2, q0 + q1 + q2);
    }
}
