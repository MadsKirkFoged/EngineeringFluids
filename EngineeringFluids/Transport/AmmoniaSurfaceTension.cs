using EngineeringUnits.Fast;
using System;

namespace EngineeringFluids.Transport;

/// <summary>
/// Liquid-vapor surface tension of ammonia, Mulero, Cachadina and Parra, J. Phys. Chem. Ref. Data 41 (2012) 043105:
/// sigma = sum a_i * (1 - T/Tc)^n_i with the correlation's own Tc = 405.4 K (CoolProp 6.4's ANCILLARIES.surface_tension).
/// </summary>
public static class AmmoniaSurfaceTension
{
    private const double Tc = 405.4; // K - the correlation's own critical temperature, not the EOS value
    private const double A0 = 0.1028, N0 = 1.211;
    private const double A1 = -0.09453, N1 = 5.585;

    /// <summary>Surface tension at a saturation temperature. NaN above the correlation's Tc (CoolProp throws there).</summary>
    public static ForcePerLength SurfaceTension(Temperature temperature)
    {
        double T = temperature.Kelvin;
        if (T > Tc)
            return ForcePerLength.FromNewtonPerMeter(double.NaN);

        double theta = 1 - T / Tc;
        return ForcePerLength.FromNewtonPerMeter(A0 * Math.Pow(theta, N0) + A1 * Math.Pow(theta, N1));
    }
}
