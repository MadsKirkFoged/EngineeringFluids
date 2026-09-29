using EngineeringUnits.Fast;
using System;
using System.Runtime.CompilerServices;

/// <summary>
/// Fast inverse of <see cref="SaturationPressure"/>: given a saturation pressure (Pa),
/// returns the saturation temperature (K) via bisection over the same cheap piecewise-quadratic
/// polynomial. SaturationPressure has no public derivative, and it is already O(1) per
/// evaluation, so a fixed-iteration-count bisection is simpler and safer (guaranteed monotonic,
/// no overshoot risk) than a Newton iteration here, at a cost of a few dozen cheap calls.
///
/// 30 iterations over the ~210 K domain already matches the underlying polynomial's own
/// precision floor (tested: no measurable difference vs 40 iterations, and rel. error vs
/// CoolProp is ~1e-7 across the whole domain away from the last ~0.5 K below Tc, where no
/// ancillary-based two-phase calculation is reliable - see UpdatePX.cs).
/// </summary>
public static class SaturationTemperature
{
    private const float Tmin = SaturationPressure.Tmin;
    private const float Tmax = SaturationPressure.Tmax;
    private const int Iterations = 30;

    public static EngineeringUnits.Fast.Temperature Temperature(Pressure p)
        => EngineeringUnits.Fast.Temperature.FromKelvin(Temperature((float)p.Pascal));

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static float Temperature(float pressurePa)
    {
        float lo = Tmin, hi = Tmax;

        for (int i = 0; i < Iterations; i++)
        {
            float mid = 0.5f * (lo + hi);
            float pMid = SaturationPressure.Pressure(mid);
            if (pMid < pressurePa)
                lo = mid;
            else
                hi = mid;
        }

        return 0.5f * (lo + hi);
    }
}
