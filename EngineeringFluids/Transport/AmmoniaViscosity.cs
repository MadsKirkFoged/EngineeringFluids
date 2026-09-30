using EngineeringUnits.Fast;
using System;

namespace EngineeringFluids.Transport;

/// <summary>
/// Viscosity of ammonia from Fenghour et al., J. Phys. Chem. Ref. Data 24 (1995) 1649, in the form and with the
/// coefficients CoolProp 6.4 uses (TRANSPORT.viscosity of its Ammonia fluid file): dilute gas + Rainwater-Friend
/// initial density dependence + modified Batschinski-Hildebrand higher-order terms. No critical enhancement.
/// </summary>
/// <remarks>
/// The correlation is empirical and written in fixed units, so it works on plain numbers internally.
/// </remarks>
public static class AmmoniaViscosity
{
    // Lennard-Jones parameters
    private const double EpsilonOverK = 386;   // K
    private const double SigmaEta = 2.957e-10; // m

    // Dilute gas, collision-integral form: eta0 = C * sqrt(M*T) / (sigma^2 * S), M in kg/kmol, sigma in nm.
    // C carries the factor of 100 missing from the published equation (see CoolProp's note on this model).
    private const double C = 2.1357e-6;
    private const double MolarMassKgPerKmol = 17.03026; // Fenghour's value - deliberately not the EOS molar mass
    private const double A0 = 4.9931822, A1 = -0.61122364, A2 = 0.0, A3 = 0.18535124, A4 = -0.11160946;

    // Rainwater-Friend second viscosity virial coefficient: B*eta = sum b_i * (T*)^(-i/2), i = 0..12
    private static readonly double[] B =
    {
        -1.7999496, 46.692621, -534.60794, 3360.4074, -13019.164, 33414.23, -58711.743,
        71426.686, -59834.012, 33652.741, -12027.35, 2434.8205, -208.07957
    };
    private const double AvogadroNumber = 6.02214129e23; // 1/mol, the value CoolProp uses here

    // Higher-order terms: sum a_i * delta^d_i * tau^t_i, delta = rho/13211.8 mol/m3, tau = 405.4 K / T
    private const double RhoReducing = 13211.8; // mol/m3
    private const double TReducing = 405.4;     // K
    private const double H0 = 4.005040600989671e-06, H1 = -1.4107915123955129e-05, H2 = 3.476074303932182e-05,
                         H3 = 4.631310990138071e-06, H4 = -3.937374461785061e-06, H5 = -1.200075068367531e-05,
                         H6 = 1.9284977991745305e-06;

    public static DynamicViscosity Viscosity(Temperature temperature, Molarity molarDensity)
    {
        double T = temperature.Kelvin;
        double rho = molarDensity.MolesPerCubicMeter;

        double tStar = T / EpsilonOverK;

        // Dilute gas
        double lnTStar = Math.Log(tStar);
        double S = Math.Exp(A0 + lnTStar * (A1 + lnTStar * (A2 + lnTStar * (A3 + lnTStar * A4))));
        double sigmaNm = SigmaEta * 1e9;
        double eta0 = C * Math.Sqrt(MolarMassKgPerKmol * T) / (sigmaNm * sigmaNm * S);

        // Initial density dependence
        double x = 1.0 / Math.Sqrt(tStar); // (T*)^(-1/2)
        double bStar = 0.0;
        for (int i = B.Length - 1; i >= 0; i--)
            bStar = bStar * x + B[i];
        double bEta = AvogadroNumber * (SigmaEta * SigmaEta * SigmaEta) * bStar; // m3/mol
        double etaInitial = eta0 * bEta * rho;

        // Higher-order terms: (d, t) = (3,0) (3,1) (2,2) (4,2) (4,3) (2,4) (4,4)
        double delta = rho / RhoReducing;
        double tau = TReducing / T;
        double d2 = delta * delta, d3 = d2 * delta, d4 = d2 * d2;
        double t2 = tau * tau, t3 = t2 * tau, t4 = t2 * t2;
        double etaHigher = H0 * d3 + H1 * d3 * tau + H2 * d2 * t2 + H3 * d4 * t2 + H4 * d4 * t3 + H5 * d2 * t4 + H6 * d4 * t4;

        return DynamicViscosity.FromPascalSecond(eta0 + etaInitial + etaHigher);
    }
}
