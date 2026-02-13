using System;

namespace EngineeringFluids;

public static class EngineeringMath
{
    /// <summary>
    /// Fast integer power for non-negative exponents.
    /// </summary>
    public static double Pow(double x, int y)
    {
        if (y < 0)
        {
            // Handle negative exponents safely
            return 1.0 / Pow(x, -y);
        }

        double result = 1.0;
        for (int i = 0; i < y; i++)
        {
            result *= x;
        }
        return result;
    }

    /// <summary>
    /// Uses integer exponentiation if exponent is (very close to) an integer.
    /// </summary>
    public static double PowWithInts(double x, double y)
    {
        // Check if y is basically an integer
        int yi = (int)Math.Round(y);
        if (Math.Abs(y - yi) < 1e-12)
        {
            return Pow(x, yi);
        }
        return Math.Pow(x, y);
    }
}
