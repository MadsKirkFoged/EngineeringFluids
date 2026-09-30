using EngineeringUnits.Fast;
using System;

namespace EngineeringFluids.Transport;

/// <summary>
/// Thermal conductivity of ammonia from Tufeu, Ivanov, Garrabos and Le Neindre, Ber. Bunsenges. Phys. Chem. 88
/// (1984) 422, as CoolProp 6.4 implements it: polynomial dilute gas + polynomial residual (TRANSPORT.conductivity of
/// its Ammonia fluid file) + the hardcoded critical enhancement (TransportRoutines::conductivity_critical_hardcoded_ammonia).
/// </summary>
/// <remarks>
/// The correlation is empirical and written in fixed units, so it works on plain numbers internally.
/// </remarks>
public static class AmmoniaThermalConductivity
{
    // Dilute gas: lambda0 = sum A_i * T^i  [W/m/K], T in K
    private const double A0 = 0.03589, A1 = -0.000175, A2 = 4.551e-07, A3 = 1.685e-10, A4 = -4.828e-13;

    // Residual: sum B_i * (rho/235)^i, i = 1..4, rho in kg/m3
    private const double B1 = 0.03808645, B2 = 0.06647986, B3 = -0.0300295, B4 = 0.00998779;
    private const double RhoReducing = 235; // kg/m3

    public static ThermalConductivity Conductivity(Temperature temperature, Density density)
    {
        double T = temperature.Kelvin;
        double rho = density.KilogramPerCubicMeter;

        double dilute = A0 + T * (A1 + T * (A2 + T * (A3 + T * A4)));

        double delta = rho / RhoReducing;
        double residual = delta * (B1 + delta * (B2 + delta * (B3 + delta * B4)));

        return ThermalConductivity.FromWattPerMeterKelvin(dilute + residual + Critical(T, rho));
    }

    // Critical enhancement, Tufeu et al. 1984 - a line-by-line port of CoolProp's version, including its constants
    // (pi to 10 digits, the 2006 Boltzmann constant). It diverges at T = 405.4 K exactly, as in CoolProp.
    private static double Critical(double T, double rho)
    {
        const double Tc = 405.4, rhoc = 235;
        const double LAMBDA = 1.2, nu = 0.63, gamma = 1.24, zeta_0_plus = 1.34e-10, a_zeta = 1, GAMMA_0_plus = 0.423e-8; // DELTA = 0.5, see below
        const double pi = 3.141592654, k_B = 1.3806504e-23;

        double t = Math.Abs((T - Tc) / Tc);
        double a_chi = a_zeta / 0.7;
        double eta_B = (2.60 + 1.6 * t) * 1e-5;
        double dPdT = (2.18 - 0.12 / Math.Exp(17.8 * t)) * 1e5; // [Pa/K]
        double X_T = 0.61 * rhoc + 16.5 * Math.Log(t);

        // Along the critical isochore (only a function of temperature) (Eq. 9)
        // CoolProp's t^DELTA with DELTA = 0.5 is a square root, and its t^-gamma / t^-nu is one power t^(nu - gamma)
        double tDelta = Math.Sqrt(t);
        double DELTA_lambda_i = LAMBDA * (k_B * T * T) / (6 * pi * eta_B * (zeta_0_plus * (1 + a_zeta * tDelta)))
                                * dPdT * dPdT * GAMMA_0_plus * Math.Pow(t, nu - gamma) * (1 + a_chi * tDelta);
        double DELTA_lambda_id = DELTA_lambda_i * Math.Exp(-36 * t * t);

        if (rho < 0.6 * rhoc)
        {
            double offset = 0.6 * rhoc - 0.96 * rhoc;
            return DELTA_lambda_id * (X_T * X_T) / (X_T * X_T + offset * offset) * (rho * rho) / ((0.6 * rhoc) * (0.6 * rhoc));
        }

        double dRho = rho - 0.96 * rhoc;
        return DELTA_lambda_id * (X_T * X_T) / (X_T * X_T + dRho * dRho);
    }
}
