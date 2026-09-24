using System;
using System.Runtime.CompilerServices;

/// <summary>
/// Fast inverse of <see cref="SaturationPressureFast"/>: given a saturation pressure (Pa),
/// returns the saturation temperature (K) via bisection over the same cheap piecewise-quadratic
/// polynomial. SaturationPressureFast has no public derivative, and it is already O(1) per
/// evaluation, so a fixed-iteration-count bisection is simpler and safer (guaranteed monotonic,
/// no overshoot risk) than a Newton iteration here, at a cost of a few dozen cheap calls.
///
/// 30 iterations over the ~210 K domain already matches the underlying polynomial's own
/// precision floor (tested: no measurable difference vs 40 iterations, and rel. error vs
/// CoolProp is ~1e-7 across the whole domain away from the last ~0.5 K below Tc, where no
/// ancillary-based two-phase calculation is reliable - see UpdatePXFast.cs).
/// </summary>
public static class SaturationTemperatureFast
{
    private const float Tmin = SaturationPressureFast.Tmin;
    private const float Tmax = SaturationPressureFast.Tmax;
    private const int Iterations = 30;

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    public static float Temperature(float pressurePa)
    {
        float lo = Tmin, hi = Tmax;

        for (int i = 0; i < Iterations; i++)
        {
            float mid = 0.5f * (lo + hi);
            float pMid = SaturationPressureFast.Pressure(mid);
            if (pMid < pressurePa)
                lo = mid;
            else
                hi = mid;
        }

        return 0.5f * (lo + hi);
    }
}
