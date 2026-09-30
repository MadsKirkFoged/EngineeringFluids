using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringFluids.Transport;
using EngineeringUnits.Fast;
using System;
using System.Threading;
using static EngineeringFluids.Helmholtz.Phase;

namespace EngineeringFluids.Fluids;

/// <remarks>
/// Every property is a function of (<see cref="Temperature"/>, <see cref="Density"/>). The Helmholtz energy and all
/// its derivatives are evaluated once, on the first property read after the state changes, and cached until
/// Temperature or Density is set again - so reading many properties costs one EOS evaluation, not one per property.
/// <para>
/// Threads: once a state is set, any number of threads may read properties of the same instance at the same time
/// (the cache is published with volatile writes, and threads that race to fill it compute identical values).
/// Changing the state - setting Temperature/Density or calling an Update method - must not overlap with anything
/// else on that instance.
/// </para>
/// </remarks>
public class Ammonia
{
    // ---------------------------
    // Primary state (single-phase)
    // ---------------------------
    private Temperature _temperature;
    private Density _density;

    public Temperature Temperature
    {
        get => _temperature;
        set { _temperature = value; InvalidateState(); }
    }

    public Density Density
    {
        get => _density;
        set { _density = value; InvalidateState(); }
    }

    // ---------------------------
    // Fixed reference points
    // ---------------------------
    // Shared by every instance (they are immutable): the solvers create short-lived Ammonia objects in their
    // inner loops, and building these per instance used to dominate the cost of `new Ammonia()`.
    public SpecificPoint Critical => CriticalPoint;
    public SpecificPoint TripleLiquid => TripleLiquidPoint;
    public SpecificPoint TripleVapor => TripleVaporPoint;
    public MolarMass MolarMass => MolarMassValue;
    public MolarEntropy GasConstant => GasConstantValue;

    private static readonly SpecificPoint CriticalPoint = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(405.56),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(21501.16668203028),
        Pressure = Pressure.FromPascal(11363400),
        MolarDensity = Molarity.FromMolesPerCubicMeter(13696),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(68.56438502935785)
    };

    private static readonly SpecificPoint TripleLiquidPoint = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(0.14111811220522047),
        Pressure = Pressure.FromPascal(6091.2231081315085),
        MolarDensity = Molarity.FromMolesPerCubicMeter(43035.33929207322),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(-1.9440525067083775e-06)
    };

    private static readonly SpecificPoint TripleVaporPoint = new SpecificPoint()
    {
        Temperature = Temperature.FromKelvin(195.495),
        MolarEnthalpy = MolarEnergy.FromJoulePerMole(25279.492873914965),
        Pressure = Pressure.FromPascal(6091.223108650368),
        MolarDensity = Molarity.FromMolesPerCubicMeter(3.763506027681136),
        MolarEntropy = MolarEntropy.FromJoulePerMoleKelvin(129.30945229032756),
    };

    private static readonly MolarMass MolarMassValue = MolarMass.FromKilogramPerMole(0.01703052);
    private static readonly MolarEntropy GasConstantValue = MolarEntropy.FromJoulePerMoleKelvin(8.3144598);

    // ---------------------------
    // Two-phase mode storage
    // ---------------------------
    private bool _isTwoPhase;
    private double _quality; // mass quality in [0,1]
    private SaturationSolver.SatResult? _satCache;
    private Ammonia? _satLiquid;
    private Ammonia? _satVapor;

    public void SetTwoPhase(SaturationSolver.SatResult sat, double q) => SetTwoPhase(sat, q, Phases.Twophase);

    // phase: the label to report - Twophase, or CriticalPoint when a flash lands exactly on it (as CoolProp does)
    internal void SetTwoPhase(SaturationSolver.SatResult sat, double q, Phases phase)
    {
        if (q < 0 || q > 1)
            throw new ArgumentOutOfRangeException(nameof(q), "Quality must be in [0,1].");

        _satCache = sat;
        _satLiquid = null;
        _satVapor = null;
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

        _phase = phase;
    }

    public void ClearTwoPhase()
    {
        _isTwoPhase = false;
        _quality = double.NaN;
        _satCache = null;
        _satLiquid = null;
        _satVapor = null;
    }

    private Ammonia SatLiquidState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");

        // A reference is published only after the object it points to is fully constructed
        Ammonia? state = Volatile.Read(ref _satLiquid);
        if (state is null)
        {
            state = new Ammonia { Temperature = _satCache.T, Density = _satCache.RhomolarL * MolarMass };
            Volatile.Write(ref _satLiquid, state);
        }
        return state;
    }

    private Ammonia SatVaporState()
    {
        if (_satCache == null)
            throw new InvalidOperationException("Two-phase state missing saturation cache.");

        Ammonia? state = Volatile.Read(ref _satVapor);
        if (state is null)
        {
            state = new Ammonia { Temperature = _satCache.T, Density = _satCache.RhomolarV * MolarMass };
            Volatile.Write(ref _satVapor, state);
        }
        return state;
    }

    // ---------------------------
    // Per-state cache
    // ---------------------------
    // Each cached value is written before its flag, and the flag with a volatile (release) write, so a thread that
    // sees the flag set also sees the value. Two threads may both compute a value; they compute the same bits.
    private readonly record struct HelmholtzState(double Tau, double Delta, IdealDerivatives Ideal, ResidualDerivatives Residual);

    private HelmholtzState _state;
    private bool _hasState;
    private double _viscositySI;
    private bool _hasViscosity;
    private double _conductivitySI;
    private bool _hasConductivity;
    private double _tsatSI;
    private bool _hasTsat;

    // The label the last flash assigned (CoolProp's per-flash conventions); null after a direct state change
    private Phases? _phase;

    private void InvalidateState()
    {
        _hasState = false;
        _hasViscosity = false;
        _hasConductivity = false;
        _hasTsat = false;
        _phase = null;
    }

    internal void SetPhase(Phases phase) => _phase = phase;

    private ref readonly HelmholtzState State
    {
        get
        {
            if (!Volatile.Read(ref _hasState))
            {
                _state = EvaluateState();
                Volatile.Write(ref _hasState, true);
            }
            return ref _state;
        }
    }

    private HelmholtzState EvaluateState()
    {
        double tau = (double)(Critical.Temperature / Temperature);
        double delta = (double)(Density / (Critical.MolarDensity * MolarMass));

        IdealDerivatives lead = IdealGasHelmholtzLead.Derivatives(delta, tau);
        IdealDerivatives logTau = IdealHelmholtzLogTau.Derivatives(delta, tau);
        IdealDerivatives planckEinstein = IdealHelmholtzPlanckEinstein.Derivatives(delta, tau);

        ResidualDerivatives power = ResidualHelmholtzPower.Derivatives(delta, tau);
        ResidualDerivatives gaussian = ResidualHelmholtzGaussian.Derivatives(delta, tau);
        ResidualDerivatives gaoB = ResidualHelmholtzGaoB.Derivatives(delta, tau);

        var ideal = new IdealDerivatives(
            lead.Value + logTau.Value + planckEinstein.Value,
            lead.dTau + logTau.dTau + planckEinstein.dTau,
            lead.dTau2 + logTau.dTau2 + planckEinstein.dTau2);

        var residual = new ResidualDerivatives(
            power.Value + gaussian.Value + gaoB.Value,
            power.dDelta + gaussian.dDelta + gaoB.dDelta,
            power.dTau + gaussian.dTau + gaoB.dTau,
            power.dDelta2 + gaussian.dDelta2 + gaoB.dDelta2,
            power.dTau2 + gaussian.dTau2 + gaoB.dTau2,
            power.dDeltadTau + gaussian.dDeltadTau + gaoB.dDeltadTau);

        return new HelmholtzState(tau, delta, ideal, residual);
    }

    // ---------------------------
    // Reduced variables
    // ---------------------------
    public double Tau => State.Tau;

    public double Delta => State.Delta;

    // ---------------------------
    // Helmholtz energy terms
    // ---------------------------
    public double Alpha0 => State.Ideal.Value;
    public double AlphaR => State.Residual.Value;
    public double Alpha => Alpha0 + AlphaR;

    public double Alpha0_dTau => State.Ideal.dTau;
    public double Alpha0_dTau2 => State.Ideal.dTau2;
    public double AlphaR_dDelta => State.Residual.dDelta;
    public double AlphaR_dTau => State.Residual.dTau;
    public double AlphaR_dDelta2 => State.Residual.dDelta2;
    public double AlphaR_dTau2 => State.Residual.dTau2;
    public double AlphaR_dDeltadTau => State.Residual.dDeltadTau;

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

    // Molar enthalpy (J/mol)
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
    // Heat capacities and speed of sound
    // ---------------------------
    // In two-phase, Cp/Cv (and the transport properties below) are evaluated at the bulk state - the
    // saturation temperature and the mixture density - exactly as CoolProp does. That point lies inside
    // the dome, so these numbers have no physical meaning there; they are kept for compatibility.

    // cv = -R * tau^2 * (alphaR_tautau + alpha0_tautau)
    public MolarEntropy MolarCv => -GasConstant * (Tau * Tau * (AlphaR_dTau2 + Alpha0_dTau2));

    // cp = cv + R * (1 + delta*alphaR_delta - delta*tau*alphaR_deltatau)^2 / (1 + 2*delta*alphaR_delta + delta^2*alphaR_deltadelta)
    public MolarEntropy MolarCp
    {
        get
        {
            double delta = Delta;
            double tau = Tau;
            double numerator = 1 + delta * AlphaR_dDelta - delta * tau * AlphaR_dDeltadTau;
            double denominator = 1 + 2 * delta * AlphaR_dDelta + delta * delta * AlphaR_dDelta2;
            return GasConstant * (-(tau * tau) * (AlphaR_dTau2 + Alpha0_dTau2) + numerator * numerator / denominator);
        }
    }

    public SpecificHeatCapacity Cv => MolarCv / MolarMass;

    public SpecificHeatCapacity Cp => MolarCp / MolarMass;

    // Like CoolProp: in two-phase only defined on the saturation lines (q = 0 or 1), where it is the saturated
    // liquid/vapor value; NaN in between, because it depends on how the phases are distributed.
    public Speed SoundSpeed
    {
        get
        {
            if (_isTwoPhase)
            {
                const double DoubleEpsilon = 2.220446049250313e-16; // C++ DBL_EPSILON, as in CoolProp
                if (Math.Abs(_quality) < DoubleEpsilon)
                    return SatLiquidState().SoundSpeed;
                if (Math.Abs(_quality - 1) < DoubleEpsilon)
                    return SatVaporState().SoundSpeed;
                return Speed.FromMeterPerSecond(double.NaN);
            }

            double delta = Delta;
            double tau = Tau;
            double numerator = 1 + delta * AlphaR_dDelta - delta * tau * AlphaR_dDeltadTau;
            double w2 = 1 + 2 * delta * AlphaR_dDelta + delta * delta * AlphaR_dDelta2
                        - numerator * numerator / (tau * tau * (AlphaR_dTau2 + Alpha0_dTau2));
            return (GasConstant * Temperature / MolarMass * w2).Sqrt();
        }
    }

    // ---------------------------
    // Transport properties
    // ---------------------------
    public DynamicViscosity DynamicViscosity
    {
        get
        {
            if (!Volatile.Read(ref _hasViscosity))
            {
                _viscositySI = AmmoniaViscosity.Viscosity(Temperature, MolarDensity).PascalSecond;
                Volatile.Write(ref _hasViscosity, true);
            }
            return DynamicViscosity.FromPascalSecond(_viscositySI);
        }
    }

    public ThermalConductivity Conductivity
    {
        get
        {
            if (!Volatile.Read(ref _hasConductivity))
            {
                _conductivitySI = AmmoniaThermalConductivity.Conductivity(Temperature, Density).WattPerMeterKelvin;
                Volatile.Write(ref _hasConductivity, true);
            }
            return ThermalConductivity.FromWattPerMeterKelvin(_conductivitySI);
        }
    }

    public double Prandtl => (double)(Cp * DynamicViscosity / Conductivity);

    // Like CoolProp: only defined in two-phase (where there is an interface); NaN otherwise.
    public ForcePerLength SurfaceTension => _isTwoPhase
        ? AmmoniaSurfaceTension.SurfaceTension(Temperature)
        : ForcePerLength.FromNewtonPerMeter(double.NaN);

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
    // Compressibility factor Z = p / (rho R T) - named as in SharpFluids
    public double Compressibility => 1 + Delta * AlphaR_dDelta;

    public MolarEnergy Gibbsmolar => GasConstant * Temperature * (1 + Alpha0 + AlphaR + Delta * AlphaR_dDelta);

    // ---------------------------
    // Phase + Quality
    // ---------------------------
    // The label follows CoolProp, whose convention depends on which update produced the state (see each Update
    // method). A state set directly through Temperature/Density gets CoolProp's plain single-phase rule.
    public Phases Phase => _phase ?? (_isTwoPhase ? Phases.Twophase : SinglePhaseLabel());

    // CoolProp's recalculate_singlephase_phase: above Pc by T vs Tc; at or below Pc, SupercriticalGas above Tc,
    // otherwise liquid or gas by density vs the critical density
    internal Phases SinglePhaseLabel()
    {
        if (Pressure > Critical.Pressure)
            return Temperature > Critical.Temperature ? Phases.Supercritical : Phases.SupercriticalLiquid;
        if (Temperature > Critical.Temperature)
            return Phases.SupercriticalGas;
        return MolarDensity > Critical.MolarDensity ? Phases.Liquid : Phases.Gas;
    }

    // Vapor mass fraction in two-phase; -1 for every single-phase state (CoolProp reports update-dependent
    // sentinels there, e.g. 10000 after a DH flash - they all mean "not two-phase").
    public double Quality => _isTwoPhase ? _quality : -1.0;

    // ---------------------------
    // Limits (the range of the Gao et al. 2020 EOS as CoolProp states it) - named as in SharpFluids
    // ---------------------------
    public Temperature CriticalTemperature => Critical.Temperature;
    public Pressure CriticalPressure => Critical.Pressure;
    public Temperature LimitTemperatureMin => TripleLiquid.Temperature;
    public Temperature LimitTemperatureMax => LimitTemperatureMaxValue;
    public Pressure LimitPressureMin => TripleLiquid.Pressure;
    public Pressure LimitPressureMax => LimitPressureMaxValue;

    // Lowest valid temperature; for a pure fluid (no brine) that is the triple point
    public Temperature T_freeze => LimitTemperatureMin;

    private static readonly Temperature LimitTemperatureMaxValue = Temperature.FromKelvin(725);
    private static readonly Pressure LimitPressureMaxValue = Pressure.FromPascal(1e9);

    // ---------------------------
    // Saturation helpers - named as in SharpFluids, solved against the EOS (not the ancillary fits)
    // ---------------------------

    /// <summary>Saturation temperature at <paramref name="pressure"/>. At or above the critical pressure this is the critical temperature.</summary>
    public Temperature GetSatTemperature(Pressure pressure)
    {
        if (pressure >= Critical.Pressure)
            return Critical.Temperature;
        if (!(pressure >= TripleLiquid.Pressure))
            throw new ArgumentOutOfRangeException(nameof(pressure), $"No saturation state below the triple-point pressure ({TripleLiquid.Pressure.Pascal} Pa). P={pressure.Pascal} Pa.");
        return this.SolveAtP(pressure).T;
    }

    /// <summary>Saturation pressure at <paramref name="temperature"/>. At or above the critical temperature this is the critical pressure.</summary>
    public Pressure GetSatPressure(Temperature temperature)
    {
        if (temperature >= Critical.Temperature)
            return Critical.Pressure;
        if (!(temperature >= TripleLiquid.Temperature))
            throw new ArgumentOutOfRangeException(nameof(temperature), $"No saturation state below the triple-point temperature ({TripleLiquid.Temperature.Kelvin} K). T={temperature.Kelvin} K.");
        return this.SaturationAtT(temperature).Psat;
    }

    /// <summary>
    /// Saturation temperature at the current pressure: the state's own temperature in two-phase, the critical
    /// temperature above the critical pressure, and NaN below the triple-point pressure.
    /// </summary>
    public Temperature Tsat
    {
        get
        {
            if (_isTwoPhase)
                return Temperature;

            if (!Volatile.Read(ref _hasTsat))
            {
                Pressure p = Pressure;
                _tsatSI = p >= TripleLiquid.Pressure ? GetSatTemperature(p).Kelvin : double.NaN;
                Volatile.Write(ref _hasTsat, true);
            }
            return Temperature.FromKelvin(_tsatSI);
        }
    }

    /// <summary>Saturated vapor density in two-phase; the state's density otherwise.</summary>
    public Density GasDensity => _isTwoPhase ? SatVaporState().Density : Density;

    /// <summary>Saturated liquid density in two-phase; the state's density otherwise.</summary>
    public Density LiquidDensity => _isTwoPhase ? SatLiquidState().Density : Density;
}
