using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Solvers;
using System;
using static EngineeringFluids.Helmholtz.Phase;

namespace EngineeringFluids.Fluids;

public class AmmoniaDouble
{
    // ---------------------------
    // Primary state (single-phase)
    // ---------------------------
    public double Temperature { get; set; }
    public double Density { get; set; }

    // ---------------------------
    // Fixed reference points
    // ---------------------------
    public readonly SpecificPointdouble Critical = new SpecificPointdouble()
    {
        Temperature = 405.56,
        MolarEnthalpy = 21501.16668203028,
        Pressure = 11363400,
        MolarDensity = 13696,
        MolarEntropy = 68.56438502935785
    };

    public readonly SpecificPointdouble TripleLiquid = new SpecificPointdouble()
    {
        Temperature = 195.495,
        MolarEnthalpy = 0.14111811220522047,
        Pressure = 6091.2231081315085,
        MolarDensity = 43035.33929207322,
        MolarEntropy = -1.9440525067083775e-06
    };

    public readonly SpecificPointdouble TripleVapor = new SpecificPointdouble()
    {
        Temperature = 195.495,
        MolarEnthalpy = 25279.492873914965,
        Pressure = 6091.223108650368,
        MolarDensity = 3.763506027681136,
        MolarEntropy = 129.30945229032756,
    };

    public readonly double MolarMass = 0.01703052;
    public readonly double GasConstant = 8.3144598;

    // ---------------------------
    // Two-phase mode storage
    // ---------------------------
    private bool _isTwoPhase;
    private double _quality; // mass quality in [0,1]
    private SaturationSolver.SatResultDouble? _satCache;

    public void SetTwoPhase(SaturationSolver.SatResultDouble sat, double q)
    {
        if (q < 0 || q > 1)
            throw new ArgumentOutOfRangeException(nameof(q), "Quality must be in [0,1].");

        _satCache = sat;
        _quality = q;
        _isTwoPhase = true;

        // Two-phase state is defined on saturation line:
        Temperature = sat.T;

        // Reporting-only mixture density via lever rule on specific volume (mass basis).
        // This density must NOT be used to evaluate single-phase EOS properties.
        double rhoL_mass = sat.RhomolarL * MolarMass; // kg/m3
        double rhoV_mass = sat.RhomolarV * MolarMass; // kg/m3

        double vMix = (1.0 - q) / rhoL_mass + q / rhoV_mass;
        Density = 1.0 / vMix;
    }

    public void ClearTwoPhase()
    {
        _isTwoPhase = false;
        _quality = double.NaN;
        _satCache = null;
    }

    private AmmoniaDouble SatLiquidState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");
        return new AmmoniaDouble
        {
            Temperature = _satCache.T,
            Density = _satCache.RhomolarL * MolarMass
        };
    }

    private AmmoniaDouble SatVaporState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");
        return new AmmoniaDouble
        {
            Temperature = _satCache.T,
            Density = _satCache.RhomolarV * MolarMass
        };
    }

    // ---------------------------
    // Reduced variables
    // ---------------------------
    public double Tau => Critical.Temperature / Temperature;

    public double Delta => Density / (Critical.MolarDensity * MolarMass);

    // ---------------------------
    // Helmholtz energy terms
    // ---------------------------
    public double Alpha0 => alpha0(Delta, Tau);
    public double AlphaR => alphaR(Delta, Tau);
    public double Alpha => Alpha0 + AlphaR;

    public double Alpha0_dTau => alpha0_dTau(Delta, Tau);
    public double AlphaR_dDelta => alphaR_dDelta(Delta, Tau);
    public double AlphaR_dTau => alphaR_dTau(Delta, Tau);
    public double AlphaR_dDelta2 => alphaR_dDelta2(Delta, Tau);

    // ---------------------------
    // Density helpers (molar)
    // ---------------------------
    public double Rhomolar => Density / MolarMass; // mol/m3

    public double MolarDensity => Density / MolarMass; // mol/m3

    // dp/drhomolar at constant T (SI: Pa / (mol/m3))
    public double dp_drhomolar_constT_SI
    {
        get
        {
            double delta = Delta;
            double term = 1.0 + 2.0 * delta * AlphaR_dDelta + (delta * delta) * AlphaR_dDelta2;
            return GasConstant * Temperature * term;
        }
    }

    // ---------------------------
    // Pressure (single-phase EOS or two-phase saturation)
    // ---------------------------
    public double Pressure
    {
        get
        {
            if (_isTwoPhase)
            {
                if (_satCache == null)
                    throw new InvalidOperationException("Two-phase state missing saturation cache.");
                return _satCache.Psat;
            }

            return MolarDensity * GasConstant * Temperature * (1 + Delta * AlphaR_dDelta);
        }
    }

    // ---------------------------
    // Thermodynamic properties
    // ---------------------------
    public double MolarEntropy => GasConstant * (Tau * (Alpha0_dTau + AlphaR_dTau) - Alpha0 - AlphaR);

    public double Entropy
    {
        get
        {
            if (_isTwoPhase)
            {
                var L = SatLiquidState();
                var V = SatVaporState();
                return (1.0 - _quality) * L.Entropy + _quality * V.Entropy;
            }

            return MolarEntropy / MolarMass;
        }
    }

    public double MolarInternalEnergy => GasConstant * Temperature * Tau * (Alpha0_dTau + AlphaR_dTau);

    public double InternalEnergy
    {
        get
        {
            if (_isTwoPhase)
            {
                var L = SatLiquidState();
                var V = SatVaporState();
                return (1.0 - _quality) * L.InternalEnergy + _quality * V.InternalEnergy;
            }

            return MolarInternalEnergy / MolarMass;
        }
    }

    // Molar enthalpy (J/mol)
    public double HMolarEnthalpy => GasConstant * Temperature * (1 + Tau * (Alpha0_dTau + AlphaR_dTau) + Delta * AlphaR_dDelta);

    public double Enthalpy
    {
        get
        {
            if (_isTwoPhase)
            {
                var L = SatLiquidState();
                var V = SatVaporState();
                return (1.0 - _quality) * L.Enthalpy + _quality * V.Enthalpy;
            }

            return HMolarEnthalpy / MolarMass;
        }
    }

    // ---------------------------
    // Fugacity (single-phase)
    // ---------------------------
    public double LNFugacityCoefficient => AlphaR + Delta * AlphaR_dDelta - Math.Log(1 + Delta * AlphaR_dDelta);

    public double FugacityCoefficient => Math.Exp(LNFugacityCoefficient);

    public double Fugacity => FugacityCoefficient * Pressure;

    // d(ln(phi))/d(delta)
    public double dLnPhi_dDelta
    {
        get
        {
            double a1 = AlphaR_dDelta;
            double a2 = AlphaR_dDelta2;
            double denom = 1.0 + Delta * a1;

            return 2.0 * a1 + Delta * a2 - (a1 + Delta * a2) / denom;
        }
    }

    // d(ln(phi))/d(rhomolar) at const T
    public double dLnPhi_dRhomolar_constT_SI => dLnPhi_dDelta / Critical.MolarDensity;

    // ---------------------------
    // Misc thermodynamic helpers (single-phase)
    // ---------------------------
    public double CompressibilityFactor => 1 + Delta * AlphaR_dDelta;

    public double Gibbsmolar => GasConstant * Temperature * (1 + Alpha0 + AlphaR + Delta * AlphaR_dDelta);

    // ---------------------------
    // Phase + Quality
    // ---------------------------
    public Phases Phase
    {
        get
        {
            if (_isTwoPhase)
                return Phases.Twophase;

            // Simple classification based on saturation pressure for subcritical
            if (Temperature < Critical.Temperature)
            {
                double psat = SaturationPressureFast.Pressure((float)Temperature);
                double p = Pressure;

                if (p > psat)
                    return Phases.Liquid;
                if (p < psat)
                    return Phases.Gas;
                return Phases.Twophase;
            }

            return Phases.Supercritical;
        }
    }

    public double Quality => _isTwoPhase ? _quality : -1.0;

    // ---------------------------
    // Helmholtz term composition
    // ---------------------------
    private static double alpha0(double delta, double tau)
    {
        return IdealGasHelmholtzLead.Alpha0(delta, tau) +
               IdealHelmholtzLogTau.Alpha0(delta, tau) +
               IdealHelmholtzPlanckEinsteinFast.Alpha0(delta, tau);
    }

    private static double alpha0_dTau(double delta, double tau)
    {
        return IdealGasHelmholtzLead.Alpha0_dTau(delta, tau) +
               IdealHelmholtzLogTau.Alpha0_dTau(delta, tau) +
               IdealHelmholtzPlanckEinsteinFast.Alpha0_dTau(delta, tau);
    }

    private static double alphaR(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR(delta, tau);
    }

    private static double alphaR_dDelta(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dDelta(delta, tau);
    }

    private static double alphaR_dTau(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dTau(delta, tau);
    }

    private static double alphaR_dDelta2(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dDelta2(delta, tau);
    }
}
