using System;

namespace EngineeringFluids.Helmholtz.Solvers;

internal static class RootFinder
{
    // Brent's method (inverse quadratic interpolation with bisection safeguard), as zbrent in Numerical Recipes.
    // Needs f(a) and f(b) of opposite sign and returns x with the bracket shrunk to xTol around it. The returned x
    // is not necessarily the last point f was evaluated at, so callers that keep state in f must evaluate f(x)
    // once more - and they check that residual themselves.
    internal static double Brent(Func<double, double> f, double a, double b, double fa, double fb, double xTol, int maxIts = 200)
    {
        if (fa == 0)
            return a;
        if (fb == 0)
            return b;
        if (!(Math.Sign(fa) != Math.Sign(fb)))
            throw new InvalidOperationException("Brent: the root is not bracketed.");

        double c = b, fc = fb, d = 0, e = 0;

        for (int iter = 0; iter < maxIts; iter++)
        {
            if (Math.Sign(fb) == Math.Sign(fc))
            {
                c = a; fc = fa;
                e = d = b - a;
            }
            if (Math.Abs(fc) < Math.Abs(fb))
            {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }

            double tol1 = 2.0 * 2.220446049250313e-16 * Math.Abs(b) + 0.5 * xTol;
            double xm = 0.5 * (c - b);
            if (Math.Abs(xm) <= tol1 || fb == 0.0)
                return b;

            if (Math.Abs(e) >= tol1 && Math.Abs(fa) > Math.Abs(fb))
            {
                double s = fb / fa, p, q;
                if (a == c)
                {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else
                {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0)
                    q = -q;
                p = Math.Abs(p);

                double min1 = 3.0 * xm * q - Math.Abs(tol1 * q);
                double min2 = Math.Abs(e * q);
                if (2.0 * p < Math.Min(min1, min2))
                {
                    e = d;
                    d = p / q;
                }
                else
                {
                    d = xm;
                    e = d;
                }
            }
            else
            {
                d = xm;
                e = d;
            }

            a = b;
            fa = fb;
            b += Math.Abs(d) > tol1 ? d : (xm >= 0 ? tol1 : -tol1);
            fb = f(b);
        }

        throw new InvalidOperationException("Brent: did not converge.");
    }
}
