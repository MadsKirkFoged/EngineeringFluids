using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringUnits;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Public API stays the same
    public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t)
        => UpdatePTCore(local, pTarget, t, Phases.Unknown, strictHint: false);

    public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint)
        => UpdatePTCore(local, pTarget, t, phaseHint, strictHint: true);


    private static void UpdatePTCore(Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint, bool strictHint)
    {
        if (local == null)
            throw new ArgumentNullException(nameof(local));
        if (pTarget == null)
            throw new ArgumentNullException(nameof(pTarget));
        if (t == null)
            throw new ArgumentNullException(nameof(t));

        // CoolProp-like PT for pure: single-phase only [1](https://coolprop.org/_static/doxygen/html/_helmholtz_e_o_s_mixture_backend_8cpp_source.html)
        local.ClearTwoPhase();
        local.Temperature = t;

        var Tc = local.Critical.Temperature;

        // Determine phase
        Phases phase;
        if (strictHint)
        {
            if (phaseHint != Phases.Liquid && phaseHint != Phases.Gas)
                throw new ArgumentException($"UpdatePT strict requires phaseHint=Liquid or Gas, got {phaseHint}");
            phase = phaseHint;
        }
        else
        {
            phase = DeterminePhaseCheapPure(local, pTarget, t);

            if (phase == Phases.Twophase)
            {
                // CoolProp behavior: PT for pure refuses two-phase [1](https://coolprop.org/_static/doxygen/html/_helmholtz_e_o_s_mixture_backend_8cpp_source.html)
                throw new InvalidOperationException(
                    $"UpdatePT: (T,P) is in/near two-phase region. Use UpdatePQ/UpdateQT or UpdatePT with phaseHint. " +
                    $"T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
            }
        }

        // Density guess (molar)
        double rhoGuessMol = GuessRhoMolar(local, pTarget, t, phase);

        // Solve molar density
        double rhoMol = SolveRhoMolar_TP_UnitState(local, pTarget, t, phase, rhoGuessMol);

        // Set final state
        local.Density = Molarity.FromMolesPerCubicMeter(rhoMol) * local.MolarMass;
    }

    private static Phases DeterminePhaseCheapPure(Ammonia a, Pressure p, Temperature t)
    {
        if (t >= a.Critical.Temperature)
            return Phases.Supercritical;

        // Fast psat(T) discrimination (no VLE solver)
        var psat = Saturation.CalculateSaturationPressure(t);

        // Tight ambiguity band around saturation line
        double rel = Math.Abs(p.Pascal - psat.Pascal) / Math.Max(psat.Pascal, 1.0);
        if (rel < 1e-6)
            return Phases.Twophase;

        return (p > psat) ? Phases.Liquid : Phases.Gas;
    }

    private static double GuessRhoMolar(Ammonia a, Pressure p, Temperature t, Phases phase)
    {
        // Ideal gas molar density (safe baseline)
        double rhoIdeal = p.Pascal / (a.GasConstant.SI * t.Kelvin);
        rhoIdeal = Math.Max(rhoIdeal, 1e-12);

        if (phase == Phases.Gas)
        {
            if (t > a.TripleLiquid.Temperature && t < a.Critical.Temperature)
            {
                var rhoV = VaporDensity.CalculateDensity(t);
                if (rhoV != null && rhoV.MolesPerCubicMeter > 0 && double.IsFinite(rhoV.MolesPerCubicMeter))
                    return rhoV.MolesPerCubicMeter;
            }
            return rhoIdeal;
        }

        if (phase == Phases.Liquid)
        {
            if (t > a.TripleLiquid.Temperature && t < a.Critical.Temperature)
            {
                var rhoL = LiquidDensity.CalculateDensity(t);
                if (rhoL != null && rhoL.MolesPerCubicMeter > 0 && double.IsFinite(rhoL.MolesPerCubicMeter))
                    return rhoL.MolesPerCubicMeter;
            }

            // fallback dense guess
            double rhoRed = a.Critical.MolarDensity.MolesPerCubicMeter;
            return Math.Max(1.5 * rhoRed, 50.0 * rhoIdeal);
        }

        // supercritical fallback
        return Math.Max(rhoIdeal, 0.5 * a.Critical.MolarDensity.MolesPerCubicMeter);
    }

    private static double SolveRhoMolar_TP_UnitState(Ammonia a, Pressure pTarget, Temperature t, Phases phase, double rhoGuessMol)
    {
        double rhoRed = a.Critical.MolarDensity.MolesPerCubicMeter;

        double rhoMin = 1e-12;
        double rhoMax = 20.0 * rhoRed;

        // Phase-aware bounds
        if (phase == Phases.Gas)
            rhoMax = Math.Min(rhoMax, 0.95 * rhoRed);
        if (phase == Phases.Liquid)
            rhoMin = Math.Max(rhoMin, 1.05 * rhoRed);

        double rho = Math.Clamp(rhoGuessMol, rhoMin, rhoMax);

        static double Resid(double pEOS, double pTarget) => (pEOS - pTarget) / pTarget;

        const double tol = 1e-10;

        // Newton first
        for (int i = 0; i < 15; i++)
        {
            a.Density = Molarity.FromMolesPerCubicMeter(rho) * a.MolarMass;

            double f = Resid(a.Pressure.Pascal, pTarget.Pascal);
            if (Math.Abs(f) < tol)
                return rho;

            double dpdrho = a.dp_drhomolar_constT_SI;
            if (!(dpdrho > 0) || !double.IsFinite(dpdrho))
                break;

            // df/drho = (1/pTarget)*dpdrho
            double step = f / (dpdrho / pTarget.Pascal);
            double rhoNew = rho - step;

            if (!double.IsFinite(rhoNew) || rhoNew <= rhoMin || rhoNew >= rhoMax)
                break;

            rho = rhoNew;
        }

        // Bisection fallback with global bracket
        double lo = rhoMin, hi = rhoMax;

        a.Density = Molarity.FromMolesPerCubicMeter(lo) * a.MolarMass;
        double fLo = Resid(a.Pressure.Pascal, pTarget.Pascal);

        a.Density = Molarity.FromMolesPerCubicMeter(hi) * a.MolarMass;
        double fHi = Resid(a.Pressure.Pascal, pTarget.Pascal);

        if (Math.Sign(fLo) == Math.Sign(fHi))
            throw new InvalidOperationException($"UpdatePT: could not bracket density root. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");

        for (int i = 0; i < 80; i++)
        {
            double mid = 0.5 * (lo + hi);

            a.Density = Molarity.FromMolesPerCubicMeter(mid) * a.MolarMass;
            double fMid = Resid(a.Pressure.Pascal, pTarget.Pascal);

            if (Math.Abs(fMid) < tol)
                return mid;

            if (Math.Sign(fMid) == Math.Sign(fLo))
            {
                lo = mid;
                fLo = fMid;
            }
            else
            {
                hi = mid;
                fHi = fMid;
            }
        }

        return 0.5 * (lo + hi);
    }

    /// <summary>
    /// Cheap phase discrimination for pure fluids using psat(T) ancillary.
    /// Returns Gas/Liquid/Supercritical/Twophase.
    /// </summary>
    private static Phases DeterminePhaseCheapPure(Ammonia a, double T, double p)
    {
        double Tc = a.Critical.Temperature.Kelvin;
        if (T >= Tc)
            return Phases.Supercritical;

        // Like CoolProp: near triple, phase determination is tricky; we still use psat(T) but keep a tight two-phase band. [1](https://coolprop.org/_static/doxygen/html/_helmholtz_e_o_s_mixture_backend_8cpp_source.html)
        double psat = Saturation.CalculateSaturationPressure(Temperature.FromKelvin(T)).Pascal;

        // Two-phase ambiguity band. Start tight; widen slightly only if you see false positives.
        double rel = Math.Abs(p - psat) / Math.Max(psat, 1.0);
        if (rel < 1e-6)
            return Phases.Twophase;

        return (p > psat) ? Phases.Liquid : Phases.Gas;
    }

    /// <summary>
    /// Phase-aware molar density guess.
    /// </summary>
    private static double GuessRhoMolar(Ammonia a, double T, double p, Phases phase)
    {
        // Ideal gas is a robust guess for vapor/supercritical-gas-like states
        double rhoIdeal = p / (a.GasConstant.SI * T);
        rhoIdeal = Math.Max(rhoIdeal, 1e-12);

        if (phase == Phases.Gas)
        {
            // Use vapor ancillary if in range (fast)
            if (T > a.TripleLiquid.Temperature.Kelvin + 1e-6 && T < a.Critical.Temperature.Kelvin - 1e-6)
            {
                var rhoV = VaporDensity.CalculateDensity(Temperature.FromKelvin(T));
                if (rhoV != null && double.IsFinite(rhoV.MolesPerCubicMeter) && rhoV.MolesPerCubicMeter > 0)
                    return Math.Max(1e-12, Math.Min(rhoV.MolesPerCubicMeter, 10.0 * a.Critical.MolarDensity.MolesPerCubicMeter));
            }
            return rhoIdeal;
        }

        if (phase == Phases.Liquid)
        {
            // Use liquid ancillary if in range (fast)
            if (T > a.TripleLiquid.Temperature.Kelvin + 1e-6 && T < a.Critical.Temperature.Kelvin - 1e-6)
            {
                var rhoL = LiquidDensity.CalculateDensity(Temperature.FromKelvin(T));
                if (rhoL != null && double.IsFinite(rhoL.MolesPerCubicMeter) && rhoL.MolesPerCubicMeter > 0)
                    return Math.Max(1e-12, Math.Min(rhoL.MolesPerCubicMeter, 20.0 * a.Critical.MolarDensity.MolesPerCubicMeter));
            }

            // Fallback: dense guess
            return Math.Max(1.5 * a.Critical.MolarDensity.MolesPerCubicMeter, 50.0 * rhoIdeal);
        }

        // Supercritical fallback
        return Math.Max(rhoIdeal, 0.5 * a.Critical.MolarDensity.MolesPerCubicMeter);
    }

    /// <summary>
    /// Solve rhomolar from (T,P) for a single-phase branch using safeguarded Newton + bisection fallback.
    /// Mirrors CoolProp's "Newton-first, bounded fallback" philosophy. [1](https://coolprop.org/_static/doxygen/html/_helmholtz_e_o_s_mixture_backend_8cpp_source.html)
    /// </summary>
    private static double SolveRhoMolar_TP(Ammonia a, double T, double pTarget, Phases phase, double rhoGuess)
    {
        double rhoRed = a.Critical.MolarDensity.MolesPerCubicMeter;

        // Global bounds
        double rhoMin = 1e-12;
        double rhoMax = 20.0 * rhoRed;

        // Phase-aware bounds to avoid crossing branches
        if (phase == Phases.Gas)
        {
            rhoMax = Math.Min(rhoMax, 0.95 * rhoRed);
        }
        else if (phase == Phases.Liquid)
        {
            rhoMin = Math.Max(rhoMin, 1.05 * rhoRed);
        }

        double rho = Math.Clamp(rhoGuess, rhoMin, rhoMax);

        // Residual scaled by pTarget (CoolProp often uses scaled residuals in solver wrappers; see DP residual example). [1](https://coolprop.org/_static/doxygen/html/_helmholtz_e_o_s_mixture_backend_8cpp_source.html)
        static double Resid(double pEOS, double pTarget) => (pEOS - pTarget) / pTarget;

        // Safeguarded Newton
        const int newtonIts = 20;
        const double tol = 1e-10;

        double lo = rhoMin, hi = rhoMax;
        double fLo = double.NaN, fHi = double.NaN;
        bool bracketReady = false;

        for (int i = 0; i < newtonIts; i++)
        {
            SetState(a, T, rho);
            double pEOS = a.Pressure.Pascal;
            double f = Resid(pEOS, pTarget);

            if (Math.Abs(f) < tol)
                return rho;

            double dpdrho = a.dp_drhomolar_constT_SI; // analytic dp/drho|T (your formula is correct)
            if (!(dpdrho > 0) || !double.IsFinite(dpdrho))
                break;

            // df/drho = (1/pTarget) * dpdrho
            double step = f / (dpdrho / pTarget);
            double rhoNew = rho - step;

            // If newton steps out of range, prepare a bracket and switch to fallback
            if (!double.IsFinite(rhoNew) || rhoNew <= rhoMin || rhoNew >= rhoMax)
            {
                bracketReady = TryBracket(a, T, pTarget, rhoMin, rhoMax, out lo, out hi, out fLo, out fHi);
                break;
            }

            // Mild damping near bounds
            if (rhoNew < rhoMin)
                rhoNew = 0.5 * (rho + rhoMin);
            if (rhoNew > rhoMax)
                rhoNew = 0.5 * (rho + rhoMax);

            rho = rhoNew;
        }

        if (!bracketReady)
        {
            bracketReady = TryBracket(a, T, pTarget, rhoMin, rhoMax, out lo, out hi, out fLo, out fHi);
        }

        if (!bracketReady)
            throw new InvalidOperationException($"UpdatePT: could not bracket density root. T={T} K, P={pTarget} Pa.");

        // Bisection fallback (simple and robust)
        for (int i = 0; i < 80; i++)
        {
            double mid = 0.5 * (lo + hi);
            SetState(a, T, mid);
            double fMid = Resid(a.Pressure.Pascal, pTarget);

            if (Math.Abs(fMid) < tol)
                return mid;

            if (Math.Sign(fMid) == Math.Sign(fLo))
            {
                lo = mid;
                fLo = fMid;
            }
            else
            {
                hi = mid;
                fHi = fMid;
            }
        }

        return 0.5 * (lo + hi);

        static void SetState(Ammonia a, double T, double rhomolar)
        {
            a.Temperature = Temperature.FromKelvin(T);
            a.Density = Molarity.FromMolesPerCubicMeter(rhomolar) * a.MolarMass;
        }

        static bool TryBracket(Ammonia a, double T, double pTarget, double rhoMin, double rhoMax,
            out double lo, out double hi, out double fLo, out double fHi)
        {
            lo = rhoMin;
            hi = rhoMax;

            SetState(a, T, lo);
            fLo = (a.Pressure.Pascal - pTarget) / pTarget;

            SetState(a, T, hi);
            fHi = (a.Pressure.Pascal - pTarget) / pTarget;

            return Math.Sign(fLo) != Math.Sign(fHi) &&
                   double.IsFinite(fLo) && double.IsFinite(fHi);
        }
    }




}
