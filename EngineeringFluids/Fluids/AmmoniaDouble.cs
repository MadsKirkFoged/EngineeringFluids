using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits;
using EngineeringUnits.Units;
using System;
using System.Collections.Concurrent;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
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

        // Optional reporting-only mixture density via lever rule on specific volume (mass basis).
        // IMPORTANT: This density must NOT be used to evaluate single-phase EOS properties.
        double rhoL_mass = sat.RhomolarL * MolarMass; // kg/m3
        double rhoV_mass = sat.RhomolarV * MolarMass; // kg/m3

        double vMix = (1.0 - q) / rhoL_mass + q / rhoV_mass;
        double rhoMix = 1.0 / vMix;
        Density = rhoMix;
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


    private void EnsureSinglePhaseState()
    {
        if (_isTwoPhase)
            return; // In two-phase mode, we don't require Density for EOS evaluation
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
    public double Tau
    {
        get
        {
            return Critical.Temperature / Temperature;
        }
    }

    public double Delta
    {
        get
        {
            EnsureSinglePhaseState(); // needs Density in single-phase
            return (Density / (Critical.MolarDensity * MolarMass));
        }
    }

    // ---------------------------
    // Helmholtz energy terms
    // ---------------------------
    public double Alpha0 => alpha0(Delta, Tau);
    public double Alpha0Fast => alpha0Fast(Delta, Tau);
    public double AlphaR => alphaR(Delta, Tau);
    public double AlphaRFast => alphaRFast(Delta, Tau);
    public double Alpha => Alpha0 + AlphaR;

    public double Alpha0_dTau => alpha0_dTau(Delta, Tau);
    public double Alpha0_dTauFast => alpha0_dTauFast(Delta, Tau);
    public double AlphaR_dDelta => alphaR_dDelta(Delta, Tau);
    public double AlphaR_dDeltaFast => alphaR_dDeltaFast(Delta, Tau);
    public double AlphaR_dTau => alphaR_dTau(Delta, Tau);
    public double AlphaR_dTauFast => alphaR_dTauFast(Delta, Tau);

    public double AlphaR_dDelta2 => alphaR_dDelta2(Delta, Tau);
    public double AlphaR_dDelta2Fast => alphaR_dDelta2Fast(Delta, Tau);

    // ---------------------------
    // Density helpers (molar)
    // ---------------------------
    public double Rhomolar
    {
        get
        {
            EnsureSinglePhaseState();
            return Density / MolarMass; // mol/m3
        }
    }

    public double MolarDensity
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
            //EnsureSinglePhaseState();

            double delta = Delta;

            double term = 1.0 + 2.0 * delta * AlphaR_dDelta + (delta * delta) * AlphaR_dDelta2;

            return GasConstant * Temperature * term;
        }
    }

    public double dp_drhomolar_constT_SIFast
    {
        get
        {
            //EnsureSinglePhaseState();

            double delta = Delta;

            double term = 1.0 + 2.0 * delta * AlphaR_dDeltaFast + (delta * delta) * AlphaR_dDelta2Fast;

            return GasConstant * Temperature * term;
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

            EnsureSinglePhaseState();
            return (MolarDensity * GasConstant * Temperature * (1 + Delta * AlphaR_dDelta));
        }
    }

    // Same formula as Pressure, but using AlphaR_dDeltaFast - needed by SolveAtTFast's
    // saturation Newton loop, where it is evaluated many times per call.
    public double PressureFast
    {
        get
        {
            if (_isTwoPhase)
            {
                if (_satCache == null)
                    throw new InvalidOperationException("Two-phase state missing saturation cache.");
                return _satCache.Psat;
            }

            EnsureSinglePhaseState();
            return (MolarDensity * GasConstant * Temperature * (1 + Delta * AlphaR_dDeltaFast));
        }
    }

    // ---------------------------
    // Thermodynamic properties
    // ---------------------------
    public double MolarEntropy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * (Tau * (Alpha0_dTau + AlphaR_dTau) - Alpha0 - AlphaR);
        }
    }

    // Same formula as MolarEntropy/Entropy, but using the Fast Alpha0/AlphaR/derivative terms
    // (see HMolarEnthalpyFast remarks) instead of the originals - needed by UpdatePSFast's
    // outer temperature search, where the plain Entropy property was the dominant per-
    // iteration cost. [benchmark-guided]
    public double MolarEntropyFast
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * (Tau * (Alpha0_dTauFast + AlphaR_dTauFast) - Alpha0Fast - AlphaRFast);
        }
    }

    public double EntropyFast => MolarEntropyFast / MolarMass;

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

    public double MolarInternalEnergy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * Temperature! * Tau * (Alpha0_dTau + AlphaR_dTau);
        }
    }

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

    // Molar enthalpy (J/mol) as "UnknownUnit" (your current approach)
    public double HMolarEnthalpy
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * Temperature! * (1 + Tau * (Alpha0_dTau + AlphaR_dTau) + Delta * AlphaR_dDelta);
        }
    }

    // Same formula as HMolarEnthalpy/Enthalpy, but using the Fast (unrolled/FMA, no per-call
    // Exp/Log savings beyond what each Fast residual class already does on its own - see
    // ResidualHelmholtz*Fast.alphaR_dTau, IdealHelmholtzPlanckEinsteinFast.Alpha0_dTau)
    // derivative terms instead of the originals. Mathematically identical, ~3x cheaper per
    // call - matters when this is evaluated repeatedly inside a solve loop (e.g. UpdatePHFast's
    // outer temperature search), where the plain Enthalpy property was the dominant per-
    // iteration cost. [benchmark-guided]
    public double HMolarEnthalpyFast
    {
        get
        {
            EnsureSinglePhaseState();
            return GasConstant * Temperature! * (1 + Tau * (Alpha0_dTauFast + AlphaR_dTauFast) + Delta * AlphaR_dDeltaFast);
        }
    }

    public double EnthalpyFast => HMolarEnthalpyFast / MolarMass;

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
    // Ancillary guesses (still useful for debugging/initialization)
    // ---------------------------
    public double MolarEnthalpyBubblePointGuess => AncillaryEnthalpyBubblePoint.MolarEnthalpy(EngineeringUnits.Temperature.FromSI(Temperature));
    public double MolarEnthalpyDewPointGuess => MolarEnthalpyBubblePointGuess + AncillaryEnthalpyDewPoint.MolarEnthalpy(EngineeringUnits.Temperature.FromSI(Temperature));

    public double EnthalpyBubblePoint => (MolarEnthalpyBubblePointGuess / MolarMass);
    public double EnthalpyDewPoint => (MolarEnthalpyDewPointGuess / MolarMass);

    public double MolarDensityBubblePoint => LiquidDensity.CalculateDensity(EngineeringUnits.Temperature.FromSI(Temperature)).SI;
    public double MolarDensityDewPoint => VaporDensity.CalculateDensity(EngineeringUnits.Temperature.FromSI(Temperature)).SI;
    public double SatPressure => Saturation.CalculateSaturationPressure(EngineeringUnits.Temperature.FromSI(Temperature)).SI;

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

    public double Fugacity => FugacityCoefficient * Pressure;

    // Same formula as LNFugacityCoefficient, but using the Fast AlphaR/AlphaR_dDelta terms -
    // needed by SolveAtTFast's saturation Newton loop.
    public double LNFugacityCoefficientFast
    {
        get
        {
            EnsureSinglePhaseState();
            return AlphaRFast + Delta * AlphaR_dDeltaFast - Math.Log(1 + Delta * AlphaR_dDeltaFast);
        }
    }

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

    public double dLnPhi_dDeltaFast
    {
        get
        {
            EnsureSinglePhaseState();

            double a1 = AlphaR_dDeltaFast;
            double a2 = AlphaR_dDelta2Fast;
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
            double rhomolar_red = Critical.MolarDensity;
            return dLnPhi_dDelta / rhomolar_red;
        }
    }

    public double dLnPhi_dRhomolar_constT_SIFast
    {
        get
        {
            EnsureSinglePhaseState();
            double rhomolar_red = Critical.MolarDensity;
            return dLnPhi_dDeltaFast / rhomolar_red;
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
            return GasConstant * Temperature * (1 + Alpha0 + AlphaR + Delta * AlphaR_dDelta);
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


            // Simple classification based on saturation pressure for subcritical
            if (Temperature < Critical.Temperature)
            {

                var psat = Saturation.CalculateSaturationPressure(EngineeringUnits.Temperature.FromSI(Temperature)).SI;
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

    private static double alpha0_dTauFast(double delta, double tau)
    {
        return IdealGasHelmholtzLead.Alpha0_dTau(delta, tau) +
               IdealHelmholtzLogTau.Alpha0_dTau(delta, tau) +
               IdealHelmholtzPlanckEinsteinFast.Alpha0_dTau(delta, tau);
    }

    private static double alpha0Fast(double delta, double tau)
    {
        return IdealGasHelmholtzLead.Alpha0(delta, tau) +
               IdealHelmholtzLogTau.Alpha0(delta, tau) +
               IdealHelmholtzPlanckEinsteinFast.Alpha0(delta, tau);
    }

    private static double alphaR(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR(delta, tau) +
               ResidualHelmholtzGaussian.alphaR(delta, tau) +
               ResidualHelmholtzGaoB.alphaR(delta, tau);
    }

    private static double alphaRFast(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR(delta, tau);
    }

    //Cache

    //private static readonly ConcurrentDictionary<(double, double), double> CachealphaR_dDelta = new();
    private static double alphaR_dDelta(double delta, double tau)
    {
            //var key = (delta, tau);

           // if (CachealphaR_dDelta.TryGetValue(key, out double local))
           //     return local;


        double result = ResidualHelmholtzPower.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dDelta(delta, tau);

          //  _ = CachealphaR_dDelta.TryAdd(key, result);

            return result;

    }
    private static double alphaR_dDeltaFast(double delta, double tau)
    {
        //var key = (delta, tau);

        // if (CachealphaR_dDelta.TryGetValue(key, out double local))
        //     return local;


        double result = ResidualHelmholtzPowerFast.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dDelta(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dDelta(delta, tau);

        //  _ = CachealphaR_dDelta.TryAdd(key, result);

        return result;

    }

    private static double alphaR_dTau(double delta, double tau)
    {
        return ResidualHelmholtzPower.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dTau(delta, tau);
    }

    private static double alphaR_dTauFast(double delta, double tau)
    {
        return ResidualHelmholtzPowerFast.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dTau(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dTau(delta, tau);
    }


    //private static readonly ConcurrentDictionary<(double, double), double> CachealphaR_dDelta2 = new();
    private static double alphaR_dDelta2(double delta, double tau)
    {
        //var key = (delta, tau);

        //if (CachealphaR_dDelta2.TryGetValue(key, out double local))
         //   return local;

        double result = ResidualHelmholtzPower.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaussian.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaoB.alphaR_dDelta2(delta, tau);


       // _ = CachealphaR_dDelta2.TryAdd(key, result);

        return result;
    }

    private static double alphaR_dDelta2Fast(double delta, double tau)
    {
        //var key = (delta, tau);

        //if (CachealphaR_dDelta2.TryGetValue(key, out double local))
        //   return local;

        double result = ResidualHelmholtzPowerFast.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaussianFast.alphaR_dDelta2(delta, tau) +
               ResidualHelmholtzGaoBFast.alphaR_dDelta2(delta, tau);


        // _ = CachealphaR_dDelta2.TryAdd(key, result);

        return result;
    }



}