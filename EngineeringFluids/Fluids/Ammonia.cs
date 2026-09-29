using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;
using static EngineeringFluids.Helmholtz.Phase;

namespace EngineeringFluids.Fluids;

public class Ammonia
{
    // ---------------------------
    // Primary state (single-phase)
    // ---------------------------
    public Temperature Temperature { get; set; }
    public Density Density { get; set; }

    // ---------------------------
    // Fixed reference points
    // ---------------------------
    public readonly SpecificPoint Critical = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(405.56),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(21501.16668203028),
        Pressure = Pressure.FromPascal(11363400),
        MolarDensity = Molarity.FromMolesPerCubicMeter(13696),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(68.56438502935785)
    };

    public readonly SpecificPoint TripleLiquid = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(0.14111811220522047),
        Pressure = Pressure.FromPascal(6091.2231081315085),
        MolarDensity = Molarity.FromMolesPerCubicMeter(43035.33929207322),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(-1.9440525067083775e-06)
    };

    public readonly SpecificPoint TripleVapor = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(25279.492873914965),
        Pressure = Pressure.FromPascal(6091.223108650368),
        MolarDensity = Molarity.FromMolesPerCubicMeter(3.763506027681136),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(129.30945229032756),
    };

    public readonly MolarMass MolarMass = MolarMass.FromKilogramPerMole(0.01703052);
    public readonly MolarEntropy GasConstant = MolarEntropy.FromJoulePerMoleKelvin(8.3144598);

    // ---------------------------
    // Two-phase mode storage
    // ---------------------------
    private bool _isTwoPhase;
    private double _quality; // mass quality in [0,1]
    private SaturationSolver.SatResult? _satCache;

    public void SetTwoPhase(SaturationSolver.SatResult sat, double q)
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
        Density rhoL_mass = sat.RhomolarL * MolarMass;
        Density rhoV_mass = sat.RhomolarV * MolarMass;

        SpecificVolume vMix = (1.0 - q) / rhoL_mass + q / rhoV_mass;
        Density = 1.0 / vMix;
    }

    public void ClearTwoPhase()
    {
        _isTwoPhase = false;
        _quality = double.NaN;
        _satCache = null;
    }

    private Ammonia SatLiquidState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");
        return new Ammonia
        {
            Temperature = _satCache.T,
            Density = _satCache.RhomolarL * MolarMass
        };
    }

    private Ammonia SatVaporState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");
        return new Ammonia
        {
            Temperature = _satCache.T,
            Density = _satCache.RhomolarV * MolarMass
        };
    }

    // ---------------------------
    // Reduced variables
    // ---------------------------
    public double Tau => (double)(Critical.Temperature / Temperature);

    public double Delta => (double)(Density / (Critical.MolarDensity * MolarMass));

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
    public Molarity Rhomolar => Density / MolarMass;

    public Molarity MolarDensity => Density / MolarMass;

    // dp/drhomolar at constant T: Pa / (mol/m3) = J/mol
    public MolarEnergy dp_drhomolar_constT
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
    public Pressure Pressure
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
    public MolarEntropy MolarEntropy => GasConstant * (Tau * (Alpha0_dTau + AlphaR_dTau) - Alpha0 - AlphaR);

    public SpecificEntropy Entropy
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

    public MolarEnergy MolarInternalEnergy => GasConstant * Temperature * Tau * (Alpha0_dTau + AlphaR_dTau);

    public SpecificEnergy InternalEnergy
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

    public MolarEnergy HMolarEnthalpy => GasConstant * Temperature * (1 + Tau * (Alpha0_dTau + AlphaR_dTau) + Delta * AlphaR_dDelta);

    public Enthalpy Enthalpy
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

    public Pressure Fugacity => FugacityCoefficient * Pressure;

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

    // d(ln(phi))/d(rhomolar) at const T: m3/mol (no named quantity, so the dimension is declared)
    [UnitDimension(BaseunitType.length, 3, BaseunitType.amountOfSubstance, -1)]
    public UnknownUnit dLnPhi_dRhomolar_constT => dLnPhi_dDelta / Critical.MolarDensity;

    // ---------------------------
    // Misc thermodynamic helpers (single-phase)
    // ---------------------------
    public double CompressibilityFactor => 1 + Delta * AlphaR_dDelta;

    public MolarEnergy Gibbsmolar => GasConstant * Temperature * (1 + Alpha0 + AlphaR + Delta * AlphaR_dDelta);

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
                Pressure psat = SaturationPressure.Pressure(Temperature);
                Pressure p = Pressure;

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
               IdealHelmholtzPlanckEinstein.Alpha0(delta, tau);
    }

    private static double alpha0_dTau(double delta, double tau)
    {
        return IdealGasHelmholtzLead.Alpha0_dTau(delta, tau) +
               IdealHelmholtzLogTau.Alpha0_dTau(delta, tau) +
               IdealHelmholtzPlanckEinstein.Alpha0_dTau(delta, tau);
    }

    private static double alphaR(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR(delta, tau) +
               ResidualHelmholtzGaussian.alphaR(delta, tau) +
               ResidualHelmholtzGaoB.alphaR(delta, tau);
    }

    private static double alphaR_dDelta(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dDelta(delta, tau);
    }

    private static double alphaR_dTau(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dTau(delta, tau);
    }

    private static double alphaR_dDelta2(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dDelta2(delta, tau);
    }
}
