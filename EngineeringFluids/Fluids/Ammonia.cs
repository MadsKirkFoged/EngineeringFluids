using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits;
using EngineeringUnits.Units;
using System;
using static EngineeringFluids.Helmholtz.Phase;

namespace EngineeringFluids.Fluids;

public class Ammonia
{


    // ---------------------------
    // Primary state (single-phase)
    // ---------------------------
    public Temperature? Temperature { get; set; }
    public Density? Density { get; set; }

    // ---------------------------
    // Fixed reference points
    // ---------------------------
    public readonly SpecificPoint Critical = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(405.56),
        MolarEnthalpy = 21501.16668203028,
        Pressure = Pressure.FromPascal(11363400),
        MolarDensity = Molarity.FromMolesPerCubicMeter(13696),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(68.56438502935785)
    };

    public readonly SpecificPoint TripleLiquid = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = 0.14111811220522047,
        Pressure = Pressure.FromPascal(6091.2231081315085),
        MolarDensity = Molarity.FromMolesPerCubicMeter(43035.33929207322),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(-1.9440525067083775e-06)
    };

    public readonly SpecificPoint TripleVapor = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = 25279.492873914965,
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

        // Optional reporting-only mixture density via lever rule on specific volume (mass basis).
        // IMPORTANT: This density must NOT be used to evaluate single-phase EOS properties.
        Density rhoL_mass = sat.RhomolarL * MolarMass; // kg/m3
        Density rhoV_mass = sat.RhomolarV * MolarMass; // kg/m3

        double vMix = (1.0 - q) / rhoL_mass.KilogramPerCubicMeter + q / rhoV_mass.KilogramPerCubicMeter;
        double rhoMix = 1.0 / vMix;
        Density = Density.FromKilogramPerCubicMeter(rhoMix);
    }

    public void ClearTwoPhase()
    {
        _isTwoPhase = false;
        _quality = double.NaN;
        _satCache = null;
    }

    // ---------------------------
    // State validation helpers
    // ---------------------------
    private void EnsureTemperature()
    {
        if (Temperature is null)
            throw new InvalidOperationException("Temperature must be set.");
    }

    private void EnsureSinglePhaseState()
    {
        EnsureTemperature();

        if (_isTwoPhase)
            return; // In two-phase mode, we don't require Density for EOS evaluation

        if (Density is null)
            throw new InvalidOperationException("Density must be set for single-phase calculations.");
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
    public double Tau
    {
        get
        {
            EnsureTemperature();
            return (double)(Critical.Temperature / Temperature);
        }
    }

    public double Delta
    {
        get
        {
            EnsureSinglePhaseState(); // needs Density in single-phase
            return (double)(Density / (Critical.MolarDensity * MolarMass));
        }
    }

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
    public Molarity Rhomolar
    {
        get
        {
            EnsureSinglePhaseState();
            return Density / MolarMass; // mol/m3
        }
    }

    public Molarity MolarDensity
    {
        get
        {
            EnsureSinglePhaseState();
            return Density / MolarMass; // mol/m3
        }
    }

    // dp/drhomolar at constant T (SI: Pa / (mol/m3))
    public double dp_drhomolar_constT_SI
    {
        get
        {
            EnsureSinglePhaseState();
            double term = 1.0
                + 2.0 * Delta * AlphaR_dDelta
                + (Delta * Delta) * AlphaR_dDelta2;

            return GasConstant.SI * Temperature!.SI * term;
        }
    }


    //public double dp_drhomolar_constT_SI
    //{
    //    get
    //    {
    //        EnsureEosCache(EosMask.AlphaR_dDelta | EosMask.AlphaR_dDelta2);

    //        double delta = _eosCache.Delta;

    //        double term = 1.0
    //            + 2.0 * delta * _eosCache.AlphaR_dDelta
    //            + (delta * delta) * _eosCache.AlphaR_dDelta2;

    //        return GasConstant.SI * Temperature!.SI * term;
    //    }
    //}


    // ---------------------------
    // Pressure (single-phase EOS or two-phase saturation)
    // ---------------------------
    public Pressure Pressure
    {
        get
        {
            EnsureTemperature();

            if (_isTwoPhase)
            {
                if (_satCache == null)
                    throw new InvalidOperationException("Two-phase state missing saturation cache.");
                return _satCache.Psat;
            }

            EnsureSinglePhaseState();
            return (MolarDensity * GasConstant * Temperature * (1 + Delta * AlphaR_dDelta));
        }
    }

    // ---------------------------
    // Thermodynamic properties
    // ---------------------------
    public MolarEntropy MolarEntropy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * (Tau * (Alpha0_dTau + AlphaR_dTau) - Alpha0 - AlphaR);
        }
    }

    public SpecificEntropy Entropy
    {
        get
        {
            EnsureTemperature();

            if (_isTwoPhase)
            {
                var L = SatLiquidState();
                var V = SatVaporState();
                return (1.0 - _quality) * L.Entropy + _quality * V.Entropy;
            }

            return MolarEntropy / MolarMass;
        }
    }

    public MolarEnergy MolarInternalEnergy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * Temperature! * Tau * (Alpha0_dTau + AlphaR_dTau);
        }
    }

    public SpecificEnergy InternalEnergy
    {
        get
        {
            EnsureTemperature();

            if (_isTwoPhase)
            {
                var L = SatLiquidState();
                var V = SatVaporState();
                return (1.0 - _quality) * L.InternalEnergy + _quality * V.InternalEnergy;
            }

            return MolarInternalEnergy / MolarMass;
        }
    }

    // Molar enthalpy (J/mol) as "UnknownUnit" (your current approach)
    public UnknownUnit HMolarEnthalpy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * Temperature! * (1 + Tau * (Alpha0_dTau + AlphaR_dTau) + Delta * AlphaR_dDelta);
        }
    }

    public Enthalpy Enthalpy
    {
        get
        {
            EnsureTemperature();

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
    // Ancillary guesses (still useful for debugging/initialization)
    // ---------------------------
    public double MolarEnthalpyBubblePointGuess => AncillaryEnthalpyBubblePoint.MolarEnthalpy(Temperature!);
    public double MolarEnthalpyDewPointGuess => MolarEnthalpyBubblePointGuess + AncillaryEnthalpyDewPoint.MolarEnthalpy(Temperature!);

    public Enthalpy EnthalpyBubblePoint => Enthalpy.FromSI(MolarEnthalpyBubblePointGuess / MolarMass.SI);
    public Enthalpy EnthalpyDewPoint => Enthalpy.FromSI(MolarEnthalpyDewPointGuess / MolarMass.SI);

    public Molarity MolarDensityBubblePoint => LiquidDensity.CalculateDensity(Temperature!);
    public Molarity MolarDensityDewPoint => VaporDensity.CalculateDensity(Temperature!);
    public Pressure SatPressure => Saturation.CalculateSaturationPressure(Temperature!);

    // ---------------------------
    // Fugacity (single-phase)
    // ---------------------------
    public double LNFugacityCoefficient
    {
        get
        {
            EnsureSinglePhaseState();
            return AlphaR + Delta * AlphaR_dDelta - Math.Log(1 + Delta * AlphaR_dDelta);
        }
    }

    public double FugacityCoefficient => Math.Exp(LNFugacityCoefficient);

    public Pressure Fugacity => FugacityCoefficient * Pressure;

    // d(ln(phi))/d(delta)
    public double dLnPhi_dDelta
    {
        get
        {
            EnsureSinglePhaseState();

            double a1 = AlphaR_dDelta;
            double a2 = AlphaR_dDelta2;
            double denom = 1.0 + Delta * a1;

            return 2.0 * a1 + Delta * a2 - (a1 + Delta * a2) / denom;
        }
    }

    // d(ln(phi))/d(rhomolar) at const T
    public double dLnPhi_dRhomolar_constT_SI
    {
        get
        {
            EnsureSinglePhaseState();
            double rhomolar_red = Critical.MolarDensity.MolesPerCubicMeter;
            return dLnPhi_dDelta / rhomolar_red;
        }
    }

    // ---------------------------
    // Misc thermodynamic helpers (single-phase)
    // ---------------------------
    public double CompressibilityFactor
    {
        get
        {
            EnsureSinglePhaseState();
            return 1 + Delta * AlphaR_dDelta;
        }
    }

    public double Gibbsmolar
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant.SI * Temperature!.SI * (1 + Alpha0 + AlphaR + Delta * AlphaR_dDelta);
        }
    }

    // ---------------------------
    // Phase + Quality
    // ---------------------------
    public Phases Phase
    {
        get
        {
            if (_isTwoPhase)
                return Phases.Twophase;

            if (Temperature is null)
                return Phases.Unknown;

            // Simple classification based on saturation pressure for subcritical
            if (Temperature < Critical.Temperature)
            {
                // If density isn't set, we can't determine state; return Unknown
                if (Density is null)
                    return Phases.Unknown;

                var psat = Saturation.CalculateSaturationPressure(Temperature);
                var p = Pressure;

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


    //// ---------------------------
    //// Cached EOS terms (single-phase only)
    //// ---------------------------
    //[Flags]
    //private enum EosMask : uint
    //{
    //    None = 0,
    //    TauDelta = 1 << 0,
    //    AlphaR_dDelta = 1 << 1,
    //    AlphaR_dDelta2 = 1 << 2,
    //    // optional later:
    //    // AlphaR        = 1 << 3,
    //    // Alpha0        = 1 << 4,
    //    // Alpha0_dTau   = 1 << 5,
    //    // AlphaR_dTau   = 1 << 6,
    //}

    //private struct EosCache
    //{
    //    public bool Valid;

    //    // Cache key
    //    public double T_K;          // numeric key only (no API change)
    //    public double Rho_kgm3;     // numeric key only (no API change)

    //    // Reduced variables
    //    public double Tau;
    //    public double Delta;

    //    // Cached derivatives (dimensionless)
    //    public double AlphaR_dDelta;
    //    public double AlphaR_dDelta2;

    //    public EosMask Mask;
    //}

    //private EosCache _eosCache;

    //private void InvalidateEosCache()
    //{
    //    _eosCache.Valid = false;
    //    _eosCache.Mask = EosMask.None;
    //}

    //private Temperature? _temperature;
    //private Density? _density;

    //public Temperature? Temperature
    //{
    //    get => _temperature;
    //    set
    //    {
    //        _temperature = value;
    //        InvalidateEosCache();
    //    }
    //}

    //public Density? Density
    //{
    //    get => _density;
    //    set
    //    {
    //        _density = value;
    //        InvalidateEosCache();
    //    }
    //}

    //private void EnsureEosCache(EosMask needed)
    //{
    //    EnsureSinglePhaseState(); // uses unit properties; OK

    //    // Cache key values (numeric only for keying, not changing API)
    //    double T_K = Temperature!.Kelvin;
    //    double rho_kgm3 = Density!.KilogramPerCubicMeter;

    //    // If state changed, reset cache
    //    if (!_eosCache.Valid || _eosCache.T_K != T_K || _eosCache.Rho_kgm3 != rho_kgm3)
    //    {
    //        _eosCache = new EosCache
    //        {
    //            Valid = true,
    //            T_K = T_K,
    //            Rho_kgm3 = rho_kgm3,
    //            Mask = EosMask.None
    //        };
    //    }

    //    // tau/delta are needed for derivatives
    //    if ((needed & EosMask.TauDelta) != 0 && (_eosCache.Mask & EosMask.TauDelta) == 0)
    //    {
    //        // Use your existing definitions to avoid changing semantics
    //        _eosCache.Tau = Tau;
    //        _eosCache.Delta = Delta;
    //        _eosCache.Mask |= EosMask.TauDelta;
    //    }

    //    // alphar_dDelta
    //    if ((needed & EosMask.AlphaR_dDelta) != 0 && (_eosCache.Mask & EosMask.AlphaR_dDelta) == 0)
    //    {
    //        EnsureEosCache(EosMask.TauDelta);

    //        double delta = _eosCache.Delta;
    //        double tau = _eosCache.Tau;

    //        _eosCache.AlphaR_dDelta =
    //            ResidualHelmholtzPower.alphaR_dDelta(delta, tau) +
    //            ResidualHelmholtzGaussian.alphaR_dDelta(delta, tau) +
    //            ResidualHelmholtzGaoB.alphaR_dDelta(delta, tau);

    //        _eosCache.Mask |= EosMask.AlphaR_dDelta;
    //    }

    //    // alphar_dDelta2
    //    if ((needed & EosMask.AlphaR_dDelta2) != 0 && (_eosCache.Mask & EosMask.AlphaR_dDelta2) == 0)
    //    {
    //        EnsureEosCache(EosMask.TauDelta);

    //        double delta = _eosCache.Delta;
    //        double tau = _eosCache.Tau;

    //        _eosCache.AlphaR_dDelta2 =
    //            ResidualHelmholtzPower.alphaR_dDelta2(delta, tau) +
    //            ResidualHelmholtzGaussian.alphaR_dDelta2(delta, tau) +
    //            ResidualHelmholtzGaoB.alphaR_dDelta2(delta, tau);

    //        _eosCache.Mask |= EosMask.AlphaR_dDelta2;
    //    }



    //}


    //public double AlphaR_dDelta
    //{
    //    get
    //    {
    //        EnsureEosCache(EosMask.AlphaR_dDelta);
    //        return _eosCache.AlphaR_dDelta;
    //    }
    //}

    //public double AlphaR_dDelta2
    //{
    //    get
    //    {
    //        EnsureEosCache(EosMask.AlphaR_dDelta2);
    //        return _eosCache.AlphaR_dDelta2;
    //    }
    //}


}