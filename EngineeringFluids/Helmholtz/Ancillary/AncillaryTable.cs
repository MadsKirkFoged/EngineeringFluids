using EngineeringUnits.Fast;
using System;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

/// <summary>
/// Shared evaluator of the saturation ancillaries <see cref="SaturationPressure"/>, <see cref="SaturationTemperature"/>,
/// <see cref="BubbleDensity"/> and <see cref="DewDensity"/>.
///
/// Each ancillary is a table of degree-6 polynomial segments over a variable x >= 0. The segment comes straight from the
/// leading bits of x's IEEE-754 representation - its exponent plus the top b mantissa bits - so there are 2^b segments per
/// octave of x, halving in width every octave toward x = 0, found with one shift (no log, no search, no branch per
/// region). With x = theta = 1 - T/Tc that puts the resolution next to the critical point, where the saturation curve
/// changes fastest. Segment 0 holds everything below the smallest octave, x in [0, 2^minExponent).
///
/// Table row per segment: [a, b, c0..c6]. f = x*a + b maps the segment onto [0, 1] and y = c0 + c1 f + ... + c6 f^6.
/// Measured (stopwatch, random T over the whole range): ~2.8 ns per call, against ~4.1 ns for the float piecewise
/// quadratics these tables replaced, at 1e-5..1e-4 better accuracy.
/// </summary>
internal static class AncillaryTable
{
    internal const int Stride = 9;

    internal const double Tc = 405.56;       // K, Gao et al. 2020
    internal const double Ttriple = 195.495; // K
    internal const double InvTc = 1 / Tc;
    internal const double ThetaMax = 1 - Ttriple / Tc;

    // theta = 1 - T/Tc, clamped to the fitted range [0, ThetaMax]; NaN stays NaN
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static double Theta(Temperature t)
    {
        double theta = Math.FusedMultiplyAdd(-t.Kelvin, InvTc, 1.0);
        return theta < 0 ? 0 : theta > ThetaMax ? ThetaMax : theta;
    }

    // Segment of x >= 0: its (exponent, top mantissa bits) key minus the key of segment 0, clamped to [0, count)
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static int Segment(double x, int shift, long firstKey, int count)
    {
        long k = (BitConverter.DoubleToInt64Bits(x) >> shift) - firstKey;
        return (ulong)k < (ulong)count ? (int)k : k < 0 ? 0 : count - 1;
    }

    // Segment must be in [0, rows) - Segment() guarantees it
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static double Polynomial(ReadOnlySpan<double> table, int segment, double x)
    {
        ref double c = ref Unsafe.Add(ref MemoryMarshal.GetReference(table), segment * Stride);
        double f = Math.FusedMultiplyAdd(x, c, Unsafe.Add(ref c, 1));
        double y = Math.FusedMultiplyAdd(Unsafe.Add(ref c, 8), f, Unsafe.Add(ref c, 7));
        y = Math.FusedMultiplyAdd(y, f, Unsafe.Add(ref c, 6));
        y = Math.FusedMultiplyAdd(y, f, Unsafe.Add(ref c, 5));
        y = Math.FusedMultiplyAdd(y, f, Unsafe.Add(ref c, 4));
        y = Math.FusedMultiplyAdd(y, f, Unsafe.Add(ref c, 3));
        return Math.FusedMultiplyAdd(y, f, Unsafe.Add(ref c, 2));
    }
}
