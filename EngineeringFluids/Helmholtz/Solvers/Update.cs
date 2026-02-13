using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static EngineeringFluids.Helmholtz.Phase;

namespace EngineeringFluids.Helmholtz.Solvers;
public static partial class Update
{


    public static void FindEnthalpyAtLiqAndGas(this Ammonia local)
    {

        Temperature TFixed = local.Temperature;

        Ammonia Liq = new Ammonia();
        Ammonia Gas = new Ammonia();

        //First guess

        Liq.Temperature = TFixed;
        Liq.Density = LiquidDensity.CalculateDensity(TFixed*0.98) * local.MolarMass;

        Gas.Temperature = TFixed;
        Gas.Density = VaporDensity.CalculateDensity(TFixed*1.02) * local.MolarMass;


        Density MaxLiq = Liq.Density;
        Density MinLiq = Gas.Density;

        Density MaxGas = Gas.Density;
        Density MinGas = Density.Zero;



        for (int i = 0; i < 100; i++)
        {
            //Liq
            Liq.Density = Liq.Density * 0.99;


            //Gas
            Gas.Density = Gas.Density * 1.01;



            //They are the same


        }


    }

    //public static void UpdatePT(this Ammonia local, Pressure p, Temperature t)
    //{

    //    local.Temperature = t;

    //    //guess a first value for the density
    //    local.Density = local.MolarDensityDewPoint * Ammonia.MolarMass;

    //    for (int i = 0; i < 100; i++)
    //    {

    //        double ratio = local.Pressure/p;

    //        //We should try with less Density
    //        local.Density /= ratio;

    //        var test2 = local.Pressure;

    //        if (ratio is < 1.00000000001 and
    //                     > 0.9999999999)
    //        {
    //            break;
    //        }

    //    }

    //}




    public static void UpdateDP(this Ammonia local,Density d, Pressure p)
    {

        local.Density = d;

        //guess a first value for the Pressure
        local.Temperature = local.Critical.Temperature;

        for (int i = 0; i < 100; i++)
        {

            double ratio = (double)(local.Pressure/p);

            //We should try with less Density
            local.Temperature /= ratio;

            var test2 = local.Pressure;

            if (ratio is < 1.00000000001 and
                         > 0.9999999999)
            {
                break;
            }

        }

    }

    public static void UpdateTQ(this Ammonia local, Temperature t, double q)
    {
        local.ClearTwoPhase();
        local.Temperature = t;

        var sat = local.SolveAtT(t);
        local.SetTwoPhase(sat, q);
    }

    public static void UpdatePQ(this Ammonia local, Pressure p, double q)
    {
        local.ClearTwoPhase();

        var sat = local.SolveAtP(p);
        local.Temperature = sat.T;
        local.SetTwoPhase(sat, q);
    }

    //public static void UpdatePH(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    //{
    //    if (local == null)
    //        throw new ArgumentNullException(nameof(local));
    //    if (pTarget == null)
    //        throw new ArgumentNullException(nameof(pTarget));
    //    if (hTarget == null)
    //        throw new ArgumentNullException(nameof(hTarget));

    //    local.ClearTwoPhase();

    //    double Pc = Ammonia.Critical.Pressure.Pascal;
    //    double Tc = Ammonia.Critical.Temperature.Kelvin;
    //    double hJ = hTarget.JoulePerKilogram;

    //    // ---------------------------------------------------------------------
    //    // A) Try saturation once (for dome detection + hinting). If it fails, continue.
    //    // ---------------------------------------------------------------------
    //    SaturationSolver.SatResult? sat = null;
    //    double TsatK = double.NaN;
    //    double hL = double.NaN, hV = double.NaN;

    //    if (pTarget.Pascal < Pc)
    //    {
    //        try
    //        {
    //            sat = SaturationSolver.SolveAtP(pTarget);
    //            TsatK = sat.T.Kelvin;

    //            var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * Ammonia.MolarMass };
    //            var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * Ammonia.MolarMass };

    //            hL = L.Enthalpy.SI;
    //            hV = V.Enthalpy.SI;

    //            double dh = hV - hL;

    //            // Near critical dh collapses; quality ill-conditioned -> treat as single-phase
    //            if (dh > 1e-3)
    //            {
    //                // Small band to avoid falling through for oracle-generated hMix
    //                double epsH = Math.Max(50.0, 1e-6 * Math.Abs(dh));
    //                if (hJ >= hL - epsH && hJ <= hV + epsH)
    //                {
    //                    double q = (hJ - hL) / dh;
    //                    q = Math.Max(0.0, Math.Min(1.0, q));
    //                    local.SetTwoPhase(sat, q);
    //                    return;
    //                }
    //            }
    //        }
    //        catch
    //        {
    //            sat = null; // continue without saturation info
    //        }
    //    }

    //    // ---------------------------------------------------------------------
    //    // B) Single-phase solve: find T such that h(P,T) = hTarget
    //    //     Do it with an imposed phase hint, but adaptively bracket only in valid region.
    //    //     This mirrors CoolProp's "flash routines + bounded solvers + phase handling". [1](https://github.com/CoolProp/CoolProp/discussions/2403)[2](https://cersonsky-lab.github.io/cbe710-notes/lecture_files/Lecture24.html)
    //    // ---------------------------------------------------------------------

    //    // Pick preferred hint if we have sat endpoints
    //    Phases preferredHint;
    //    if (sat != null)
    //    {
    //        if (hJ < hL)
    //            preferredHint = Phases.Liquid;
    //        else if (hJ > hV)
    //            preferredHint = Phases.Gas;
    //        else
    //            preferredHint = Phases.Unknown; // should not happen (handled above), but safe
    //    }
    //    else
    //    {
    //        preferredHint = Phases.Gas; // good default; we'll retry others
    //    }

    //    Phases[] tries;
    //    if (preferredHint == Phases.Gas)
    //        tries = new[] { Phases.Gas, Phases.Liquid, Phases.Unknown };
    //    else if (preferredHint == Phases.Liquid)
    //        tries = new[] { Phases.Liquid, Phases.Gas, Phases.Unknown };
    //    else
    //        tries = new[] { Phases.Unknown, Phases.Gas, Phases.Liquid };

    //    double TminGlobal = Math.Max(Ammonia.TripleLiquid.Temperature.Kelvin + 1e-3, 150.0);
    //    double TmaxGlobal = 2500.0;

    //    // TryEval: safe evaluation of f(T)=h(P,T)-hTarget for a given hint
    //    bool TryEval(double TK, Phases hint, out double f)
    //    {
    //        f = double.NaN;
    //        try
    //        {
    //            // If we know Tsat, keep T on correct side for liquid/vapor
    //            if (sat != null)
    //            {
    //                if (hint == Phases.Liquid)
    //                    TK = Math.Min(TK, TsatK - 1e-6);
    //                if (hint == Phases.Gas)
    //                    TK = Math.Max(TK, TsatK + 1e-6);
    //            }

    //            local.UpdatePT(pTarget, Temperature.FromKelvin(TK), hint);
    //            f = local.Enthalpy.SI - hJ;
    //            return double.IsFinite(f);
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }

    //    // Find a valid starting temperature for a hint (where UpdatePT succeeds)
    //    bool FindValidStart(Phases hint, out double Tstart, out double fstart)
    //    {
    //        Tstart = double.NaN;
    //        fstart = double.NaN;

    //        // Starting guess:
    //        // - If Tsat is known, start just on the correct side
    //        // - Else start near Tc (works well broadly)
    //        double T0;
    //        if (sat != null)
    //        {
    //            if (hint == Phases.Gas)
    //                T0 = TsatK + 5.0;
    //            else if (hint == Phases.Liquid)
    //                T0 = TsatK - 5.0;
    //            else
    //                T0 = TsatK + 10.0;
    //        }
    //        else
    //        {
    //            T0 = local.Temperature?.Kelvin ?? (Tc + 20.0);
    //        }

    //        T0 = Math.Max(TminGlobal, Math.Min(TmaxGlobal, T0));

    //        // Try T0, then walk outward depending on hint
    //        double T = T0;
    //        double step = 10.0;

    //        for (int k = 0; k < 60; k++)
    //        {
    //            if (TryEval(T, hint, out double f))
    //            {
    //                Tstart = T;
    //                fstart = f;
    //                return true;
    //            }

    //            if (hint == Phases.Gas)
    //            {
    //                // Move upward until gas solve becomes feasible
    //                T = Math.Min(TmaxGlobal, T + step);
    //            }
    //            else if (hint == Phases.Liquid)
    //            {
    //                // Move downward until liquid solve becomes feasible
    //                T = Math.Max(TminGlobal, T - step);
    //            }
    //            else
    //            {
    //                // Unknown: explore both sides by increasing T
    //                T = Math.Min(TmaxGlobal, T + step);
    //            }

    //            step *= 1.2;

    //            if (T <= TminGlobal + 1e-9 || T >= TmaxGlobal - 1e-9)
    //                break;
    //        }

    //        return false;
    //    }

    //    // Bracket and solve for one hint
    //    bool TrySolveForHint(Phases hint)
    //    {
    //        if (!FindValidStart(hint, out double T0, out double f0))
    //            return false;

    //        // If already close enough
    //        if (Math.Abs(f0) < 1e-3)
    //        {
    //            local.ClearTwoPhase();
    //            return true;
    //        }

    //        // Now expand to bracket sign change
    //        double Tlo = T0, Thi = T0;
    //        double flo = f0, fhi = f0;

    //        double step = 20.0;

    //        for (int k = 0; k < 120 && Math.Sign(flo) == Math.Sign(fhi); k++)
    //        {
    //            // Expand upward and downward, preferring direction based on sign
    //            // For monotone h(T) in a given phase, this will bracket quickly.
    //            double Ttry;

    //            if (f0 < 0)
    //            {
    //                // Need higher T to raise h
    //                Ttry = Math.Min(TmaxGlobal, Thi + step);
    //                if (TryEval(Ttry, hint, out double ftry))
    //                {
    //                    Thi = Ttry;
    //                    fhi = ftry;
    //                }
    //                else
    //                {
    //                    // If we can't evaluate upward, give up this hint
    //                    break;
    //                }
    //            }
    //            else
    //            {
    //                // Need lower T to reduce h
    //                Ttry = Math.Max(TminGlobal, Tlo - step);
    //                if (TryEval(Ttry, hint, out double ftry))
    //                {
    //                    Tlo = Ttry;
    //                    flo = ftry;
    //                }
    //                else
    //                {
    //                    // If we can't evaluate downward, give up this hint
    //                    break;
    //                }
    //            }

    //            step *= 1.25;
    //        }

    //        if (Math.Sign(flo) == Math.Sign(fhi))
    //            return false;

    //        // Bisection inside valid bracket
    //        const int maxIts = 220;
    //        const double absTolH = 1e-3;
    //        const double relTolH = 1e-10;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double Tm = 0.5 * (Tlo + Thi);

    //            if (!TryEval(Tm, hint, out double fm))
    //            {
    //                // If midpoint fails, nudge towards the side that evaluates
    //                double Tm1 = 0.5 * (Tlo + Tm);
    //                double Tm2 = 0.5 * (Tm + Thi);

    //                bool ok1 = TryEval(Tm1, hint, out double f1);
    //                bool ok2 = TryEval(Tm2, hint, out double f2);

    //                if (ok1 && !ok2)
    //                { Tm = Tm1; fm = f1; }
    //                else if (!ok1 && ok2)
    //                { Tm = Tm2; fm = f2; }
    //                else if (ok1 && ok2)
    //                {
    //                    // pick the one closer to zero
    //                    if (Math.Abs(f1) < Math.Abs(f2))
    //                    { Tm = Tm1; fm = f1; }
    //                    else
    //                    { Tm = Tm2; fm = f2; }
    //                }
    //                else
    //                {
    //                    return false;
    //                }
    //            }

    //            if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
    //            {
    //                // Finalize state at Tm
    //                local.UpdatePT(pTarget, Temperature.FromKelvin(Tm), hint);
    //                local.ClearTwoPhase();
    //                return true;
    //            }

    //            if (Math.Sign(fm) == Math.Sign(flo))
    //            {
    //                Tlo = Tm;
    //                flo = fm;
    //            }
    //            else
    //            {
    //                Thi = Tm;
    //                fhi = fm;
    //            }
    //        }

    //        return false;
    //    }

    //    foreach (var hint in tries)
    //    {
    //        if (TrySolveForHint(hint))
    //            return;
    //    }

    //    throw new InvalidOperationException(
    //        $"UpdatePH failed. Could not solve single-phase at P={pTarget.Pascal} Pa, h={hJ} J/kg. " +
    //        $"Tried hints: {string.Join(",", tries)} " +
    //        (sat != null ? $"(Tsat={TsatK}K, hL={hL}, hV={hV})" : "(no saturation info)"));
    //}


    //public static void UpdatePH(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    //{
    //    if (local == null)
    //        throw new ArgumentNullException(nameof(local));
    //    if (pTarget == null)
    //        throw new ArgumentNullException(nameof(pTarget));
    //    if (hTarget == null)
    //        throw new ArgumentNullException(nameof(hTarget));

    //    local.ClearTwoPhase();

    //    double Pc = Ammonia.Critical.Pressure.Pascal;
    //    double Tc = Ammonia.Critical.Temperature.Kelvin;
    //    double hJ = hTarget.JoulePerKilogram;

    //    // ---------------------------------------------------------------------
    //    // A) Try saturation once (subcritical only) for dome detection + hinting
    //    // ---------------------------------------------------------------------
    //    SaturationSolver.SatResult? sat = null;
    //    double TsatK = double.NaN;
    //    double hL = double.NaN, hV = double.NaN;

    //    if (pTarget.Pascal < Pc)
    //    {
    //        try
    //        {
    //            sat = SaturationSolver.SolveAtP(pTarget);
    //            TsatK = sat.T.Kelvin;

    //            var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * Ammonia.MolarMass };
    //            var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * Ammonia.MolarMass };

    //            hL = L.Enthalpy.SI;
    //            hV = V.Enthalpy.SI;

    //            double dh = hV - hL;

    //            // Near critical dh collapses; quality ill-conditioned -> treat as single-phase
    //            if (dh > 1e-3)
    //            {
    //                // Buffer so CoolProp-based hMix doesn't "just miss"
    //                double epsH = Math.Max(50.0, 1e-6 * Math.Abs(dh));

    //                if (hJ >= hL - epsH && hJ <= hV + epsH)
    //                {
    //                    double q = (hJ - hL) / dh;
    //                    q = Math.Max(0.0, Math.Min(1.0, q));
    //                    local.SetTwoPhase(sat, q);
    //                    return;
    //                }
    //            }
    //        }
    //        catch
    //        {
    //            sat = null; // fall back to blind single-phase solve
    //        }
    //    }

    //    // ---------------------------------------------------------------------
    //    // B) Single-phase solve helpers
    //    // ---------------------------------------------------------------------

    //    // 1) Residual using your robust UNHINTED UpdatePT(P,T)
    //    bool TryEvalBlind(double TK, out double f)
    //    {
    //        f = double.NaN;
    //        try
    //        {
    //            local.UpdatePT(pTarget, Temperature.FromKelvin(TK));
    //            f = local.Enthalpy.SI - hJ;
    //            return double.IsFinite(f);
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }

    //    // 2) Residual using hinted UpdatePT(P,T,phaseHint) (strict)
    //    bool TryEvalHint(double TK, Phases hint, out double f)
    //    {
    //        f = double.NaN;
    //        try
    //        {
    //            if (sat != null)
    //            {
    //                if (hint == Phases.Liquid)
    //                    TK = Math.Min(TK, TsatK - 1e-6);
    //                if (hint == Phases.Gas)
    //                    TK = Math.Max(TK, TsatK + 1e-6);
    //            }

    //            local.UpdatePT(pTarget, Temperature.FromKelvin(TK), hint);
    //            f = local.Enthalpy.SI - hJ;
    //            return double.IsFinite(f);
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }

    //    double Tmin = Math.Max(Ammonia.TripleLiquid.Temperature.Kelvin + 1e-3, 150.0);
    //    double Tmax = 2500.0;

    //    // ---------------------------------------------------------------------
    //    // C) First attempt: BLIND single-phase root solve (no saturation required)
    //    // ---------------------------------------------------------------------
    //    bool TrySolveBlind(out double Tsol)
    //    {
    //        Tsol = double.NaN;

    //        // Choose an initial guess near previous T or near Tc
    //        double T0 = local.Temperature?.Kelvin ?? (Tc + 20.0);
    //        T0 = Math.Max(Tmin, Math.Min(Tmax, T0));

    //        if (!TryEvalBlind(T0, out double f0))
    //            f0 = double.NaN;

    //        // Find a valid start by scanning upward then downward
    //        if (!double.IsFinite(f0))
    //        {
    //            double step = 20.0;
    //            bool found = false;

    //            for (int k = 0; k < 60; k++)
    //            {
    //                double Tup = Math.Min(Tmax, T0 + step);
    //                if (TryEvalBlind(Tup, out f0))
    //                {
    //                    T0 = Tup;
    //                    found = true;
    //                    break;
    //                }

    //                double Tdn = Math.Max(Tmin, T0 - step);
    //                if (TryEvalBlind(Tdn, out f0))
    //                {
    //                    T0 = Tdn;
    //                    found = true;
    //                    break;
    //                }

    //                step *= 1.25;
    //                if (Tup >= Tmax - 1e-9 && Tdn <= Tmin + 1e-9)
    //                    break;
    //            }

    //            if (!found)
    //                return false;
    //        }

    //        // Now bracket: expand both directions until sign change
    //        double Tlo = T0, Thi = T0;
    //        double flo = f0, fhi = f0;

    //        double step2 = 25.0;
    //        for (int k = 0; k < 140 && Math.Sign(flo) == Math.Sign(fhi); k++)
    //        {
    //            bool moved = false;

    //            double Tup = Math.Min(Tmax, Thi + step2);
    //            if (Tup > Thi + 1e-12 && TryEvalBlind(Tup, out double fup))
    //            {
    //                Thi = Tup;
    //                fhi = fup;
    //                moved = true;
    //            }

    //            double Tdn = Math.Max(Tmin, Tlo - step2);
    //            if (Tdn < Tlo - 1e-12 && TryEvalBlind(Tdn, out double fdn))
    //            {
    //                Tlo = Tdn;
    //                flo = fdn;
    //                moved = true;
    //            }

    //            if (!moved)
    //                break;
    //            step2 *= 1.18;
    //        }

    //        if (Math.Sign(flo) == Math.Sign(fhi))
    //            return false;

    //        // Bisection
    //        const int maxIts = 240;
    //        const double absTolH = 1e-3;
    //        const double relTolH = 1e-10;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double Tm = 0.5 * (Tlo + Thi);
    //            if (!TryEvalBlind(Tm, out double fm))
    //                return false;

    //            if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
    //            {
    //                Tsol = Tm;
    //                return true;
    //            }

    //            if (Math.Sign(fm) == Math.Sign(flo))
    //            {
    //                Tlo = Tm;
    //                flo = fm;
    //            }
    //            else
    //            {
    //                Thi = Tm;
    //                fhi = fm;
    //            }
    //        }

    //        return false;
    //    }

    //    if (TrySolveBlind(out double Tblind))
    //    {
    //        local.UpdatePT(pTarget, Temperature.FromKelvin(Tblind));
    //        local.ClearTwoPhase();
    //        return;
    //    }

    //    // ---------------------------------------------------------------------
    //    // D) Second attempt: hinted solves (only helpful if saturation exists)
    //    //     (CoolProp uses different routines/strategies; here we try hints explicitly) [1](https://cersonsky-lab.github.io/cbe710-notes/lecture_files/Lecture24.html)[2](https://github.com/CoolProp/CoolProp/discussions/2403)
    //    // ---------------------------------------------------------------------
    //    Phases preferredHint = Phases.Unknown;

    //    if (sat != null)
    //    {
    //        if (hJ < hL)
    //            preferredHint = Phases.Liquid;
    //        else if (hJ > hV)
    //            preferredHint = Phases.Gas;
    //        else
    //            preferredHint = Phases.Unknown;
    //    }
    //    else
    //    {
    //        preferredHint = Phases.Unknown;
    //    }

    //    Phases[] tries;
    //    if (preferredHint == Phases.Gas)
    //        tries = new[] { Phases.Gas, Phases.Liquid, Phases.Unknown };
    //    else if (preferredHint == Phases.Liquid)
    //        tries = new[] { Phases.Liquid, Phases.Gas, Phases.Unknown };
    //    else
    //        tries = new[] { Phases.Unknown, Phases.Gas, Phases.Liquid };

    //    bool TrySolveHint(Phases hint, out double Tsol2)
    //    {
    //        Tsol2 = double.NaN;

    //        double TminWork = Tmin;
    //        double TmaxWork = Tmax;

    //        if (sat != null)
    //        {
    //            if (hint == Phases.Liquid)
    //                TmaxWork = Math.Min(TmaxWork, TsatK - 1e-6);
    //            if (hint == Phases.Gas)
    //                TminWork = Math.Max(TminWork, TsatK + 1e-6);
    //        }

    //        if (TmaxWork <= TminWork + 1e-9)
    //            return false;

    //        if (!TryEvalHint(TminWork, hint, out double fmin))
    //            return false;
    //        if (!TryEvalHint(TmaxWork, hint, out double fmax))
    //            return false;

    //        int expand = 0;
    //        while (Math.Sign(fmin) == Math.Sign(fmax) && expand < 12)
    //        {
    //            if (hint == Phases.Gas)
    //            {
    //                TmaxWork *= 1.5;
    //                if (TmaxWork > 3000.0)
    //                    break;
    //                if (!TryEvalHint(TmaxWork, hint, out fmax))
    //                    return false;
    //            }
    //            else if (hint == Phases.Liquid)
    //            {
    //                TminWork = Math.Max(150.0, TminWork - 25.0);
    //                if (!TryEvalHint(TminWork, hint, out fmin))
    //                    return false;
    //            }
    //            else
    //            {
    //                TmaxWork *= 1.5;
    //                if (TmaxWork > 3000.0)
    //                    break;
    //                if (!TryEvalHint(TmaxWork, hint, out fmax))
    //                    return false;
    //            }
    //            expand++;
    //        }

    //        if (Math.Sign(fmin) == Math.Sign(fmax))
    //            return false;

    //        const int maxIts = 220;
    //        const double absTolH = 1e-3;
    //        const double relTolH = 1e-10;

    //        double a = TminWork, b = TmaxWork;
    //        double fa = fmin, fb = fmax;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double m = 0.5 * (a + b);
    //            if (!TryEvalHint(m, hint, out double fm))
    //                return false;

    //            if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
    //            {
    //                Tsol2 = m;
    //                return true;
    //            }

    //            if (Math.Sign(fm) == Math.Sign(fa))
    //            {
    //                a = m;
    //                fa = fm;
    //            }
    //            else
    //            {
    //                b = m;
    //                fb = fm;
    //            }
    //        }

    //        return false;
    //    }

    //    foreach (var hint in tries)
    //    {
    //        if (TrySolveHint(hint, out double Tsol3))
    //        {
    //            if (hint == Phases.Unknown)
    //                local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol3));
    //            else
    //                local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol3), hint);

    //            local.ClearTwoPhase();
    //            return;
    //        }
    //    }

    //    throw new InvalidOperationException(
    //        $"UpdatePH failed. Could not solve single-phase at P={pTarget.Pascal} Pa, h={hJ} J/kg. " +
    //        $"Tried hints: {string.Join(",", tries)} " +
    //        (sat != null ? $"(Tsat={TsatK}K, hL={hL}, hV={hV})" : "(no saturation info)"));
    //}

    //public static void UpdatePH(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    //{
    //    if (local == null)
    //        throw new ArgumentNullException(nameof(local));
    //    if (pTarget == null)
    //        throw new ArgumentNullException(nameof(pTarget));
    //    if (hTarget == null)
    //        throw new ArgumentNullException(nameof(hTarget));

    //    local.ClearTwoPhase();

    //    double Pc = Ammonia.Critical.Pressure.Pascal;
    //    double Tc = Ammonia.Critical.Temperature.Kelvin;
    //    double hJ = hTarget.JoulePerKilogram;

    //    // ---------------------------------------------------------------------
    //    // A) Two-phase detection (subcritical only)
    //    //    Strategy:
    //    //      1) Try sat solve at P (SolveAtP)
    //    //      2) If that fails, compute Tsat(P) from ancillary and try SolveAtT(Tsat)
    //    //      3) If saturation exists, compute hL/hV from your EOS and check hJ in [hL,hV] (+ buffer)
    //    // ---------------------------------------------------------------------
    //    SaturationSolver.SatResult? sat = null;
    //    double TsatK = double.NaN;
    //    double hL = double.NaN, hV = double.NaN;

    //    if (pTarget.Pascal < Pc)
    //    {
    //        // Try saturation-at-pressure
    //        try
    //        {
    //            sat = SaturationSolver.SolveAtP(pTarget);
    //        }
    //        catch
    //        {
    //            sat = null;
    //        }

    //        // Fallback: Tsat(P) from ancillary, then sat-at-temperature if available
    //        if (sat == null)
    //        {
    //            try
    //            {
    //                var Tsat = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
    //                TsatK = Tsat.Kelvin;

    //                // Fallback attempt: compute saturation state at Tsat
    //                // (many solvers are more stable at fixed T than fixed P)
    //                sat = SaturationSolver.SolveAtT(Tsat);
    //            }
    //            catch
    //            {
    //                sat = null;
    //            }
    //        }

    //        if (sat != null)
    //        {
    //            TsatK = sat.T.Kelvin;

    //            try
    //            {
    //                // Build saturated endpoints using your EOS at (Tsat, rho_sat)
    //                var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * Ammonia.MolarMass };
    //                var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * Ammonia.MolarMass };

    //                hL = L.Enthalpy.SI;
    //                hV = V.Enthalpy.SI;

    //                double dh = hV - hL;

    //                // If latent heat collapses (near critical), quality becomes ill-conditioned -> treat as single-phase
    //                if (dh > 1e-3)
    //                {
    //                    // Buffer to avoid falling through due to small differences between your EOS endpoints and CoolProp endpoints
    //                    double epsH = Math.Max(50.0, 1e-6 * Math.Abs(dh));

    //                    if (hJ >= hL - epsH && hJ <= hV + epsH)
    //                    {
    //                        double q = (hJ - hL) / dh;
    //                        q = Math.Max(0.0, Math.Min(1.0, q));
    //                        local.SetTwoPhase(sat, q);
    //                        return;
    //                    }
    //                }
    //            }
    //            catch
    //            {
    //                // If we cannot compute hL/hV, fall through to single-phase.
    //                sat = null;
    //            }
    //        }
    //    }

    //    // ---------------------------------------------------------------------
    //    // B) Single-phase solve utilities
    //    // ---------------------------------------------------------------------

    //    double Tmin = Math.Max(Ammonia.TripleLiquid.Temperature.Kelvin + 1e-3, 150.0);
    //    double Tmax = 2500.0;

    //    // Residual using robust unhinted UpdatePT(P,T)
    //    bool TryEvalBlind(double TK, out double f)
    //    {
    //        f = double.NaN;
    //        try
    //        {
    //            local.UpdatePT(pTarget, Temperature.FromKelvin(TK));
    //            f = local.Enthalpy.SI - hJ;
    //            return double.IsFinite(f);
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }

    //    // Residual using strict hinted UpdatePT(P,T,hint)
    //    bool TryEvalHint(double TK, Phases hint, out double f)
    //    {
    //        f = double.NaN;
    //        try
    //        {
    //            // If saturation known, keep on the proper side to avoid crossing into two-phase
    //            if (sat != null)
    //            {
    //                if (hint == Phases.Liquid)
    //                    TK = Math.Min(TK, TsatK - 1e-6);
    //                if (hint == Phases.Gas)
    //                    TK = Math.Max(TK, TsatK + 1e-6);
    //            }

    //            local.UpdatePT(pTarget, Temperature.FromKelvin(TK), hint);
    //            f = local.Enthalpy.SI - hJ;
    //            return double.IsFinite(f);
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }

    //    // ---------------------------------------------------------------------
    //    // C) First attempt: BLIND single-phase flash (no saturation required)
    //    // ---------------------------------------------------------------------
    //    bool TrySolveBlind(out double Tsol)
    //    {
    //        Tsol = double.NaN;

    //        // Initial guess near Tc or last known temperature
    //        double T0 = local.Temperature?.Kelvin ?? (Tc + 20.0);
    //        T0 = Math.Max(Tmin, Math.Min(Tmax, T0));

    //        if (!TryEvalBlind(T0, out double f0))
    //        {
    //            // scan for a valid start
    //            double step = 20.0;
    //            bool found = false;
    //            for (int k = 0; k < 60; k++)
    //            {
    //                double Tup = Math.Min(Tmax, T0 + step);
    //                if (TryEvalBlind(Tup, out f0))
    //                { T0 = Tup; found = true; break; }

    //                double Tdn = Math.Max(Tmin, T0 - step);
    //                if (TryEvalBlind(Tdn, out f0))
    //                { T0 = Tdn; found = true; break; }

    //                step *= 1.25;
    //                if (Tup >= Tmax - 1e-9 && Tdn <= Tmin + 1e-9)
    //                    break;
    //            }
    //            if (!found)
    //                return false;
    //        }

    //        // bracket by expanding both directions
    //        double Tlo = T0, Thi = T0;
    //        double flo = f0, fhi = f0;

    //        double step2 = 25.0;
    //        for (int k = 0; k < 140 && Math.Sign(flo) == Math.Sign(fhi); k++)
    //        {
    //            bool moved = false;

    //            double Tup = Math.Min(Tmax, Thi + step2);
    //            if (Tup > Thi + 1e-12 && TryEvalBlind(Tup, out double fup))
    //            {
    //                Thi = Tup;
    //                fhi = fup;
    //                moved = true;
    //            }

    //            double Tdn = Math.Max(Tmin, Tlo - step2);
    //            if (Tdn < Tlo - 1e-12 && TryEvalBlind(Tdn, out double fdn))
    //            {
    //                Tlo = Tdn;
    //                flo = fdn;
    //                moved = true;
    //            }

    //            if (!moved)
    //                break;
    //            step2 *= 1.18;
    //        }

    //        if (Math.Sign(flo) == Math.Sign(fhi))
    //            return false;

    //        // bisection
    //        const int maxIts = 240;
    //        const double absTolH = 1e-3;
    //        const double relTolH = 1e-10;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double Tm = 0.5 * (Tlo + Thi);
    //            if (!TryEvalBlind(Tm, out double fm))
    //                return false;

    //            if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
    //            {
    //                Tsol = Tm;
    //                return true;
    //            }

    //            if (Math.Sign(fm) == Math.Sign(flo))
    //            {
    //                Tlo = Tm;
    //                flo = fm;
    //            }
    //            else
    //            {
    //                Thi = Tm;
    //                fhi = fm;
    //            }
    //        }

    //        return false;
    //    }

    //    if (TrySolveBlind(out double Tblind))
    //    {
    //        local.UpdatePT(pTarget, Temperature.FromKelvin(Tblind));
    //        local.ClearTwoPhase();
    //        return;
    //    }

    //    // ---------------------------------------------------------------------
    //    // D) Second attempt: hinted flash (CoolProp-like "try other strategy") [2](https://cersonsky-lab.github.io/cbe710-notes/lecture_files/Lecture24.html)[1](https://github.com/CoolProp/CoolProp/discussions/2403)
    //    // ---------------------------------------------------------------------
    //    Phases preferredHint = Phases.Unknown;

    //    if (sat != null)
    //    {
    //        if (hJ < hL)
    //            preferredHint = Phases.Liquid;
    //        else if (hJ > hV)
    //            preferredHint = Phases.Gas;
    //        else
    //            preferredHint = Phases.Unknown;
    //    }
    //    else
    //    {
    //        preferredHint = Phases.Unknown;
    //    }

    //    Phases[] tries;
    //    if (preferredHint == Phases.Gas)
    //        tries = new[] { Phases.Gas, Phases.Liquid, Phases.Unknown };
    //    else if (preferredHint == Phases.Liquid)
    //        tries = new[] { Phases.Liquid, Phases.Gas, Phases.Unknown };
    //    else
    //        tries = new[] { Phases.Unknown, Phases.Gas, Phases.Liquid };

    //    bool TrySolveHint(Phases hint, out double Tsol2)
    //    {
    //        Tsol2 = double.NaN;

    //        double TminWork = Tmin;
    //        double TmaxWork = Tmax;

    //        if (sat != null)
    //        {
    //            if (hint == Phases.Liquid)
    //                TmaxWork = Math.Min(TmaxWork, TsatK - 1e-6);
    //            if (hint == Phases.Gas)
    //                TminWork = Math.Max(TminWork, TsatK + 1e-6);
    //        }

    //        if (TmaxWork <= TminWork + 1e-9)
    //            return false;

    //        if (!TryEvalHint(TminWork, hint, out double fmin))
    //            return false;
    //        if (!TryEvalHint(TmaxWork, hint, out double fmax))
    //            return false;

    //        int expand = 0;
    //        while (Math.Sign(fmin) == Math.Sign(fmax) && expand < 12)
    //        {
    //            if (hint == Phases.Gas)
    //            {
    //                TmaxWork *= 1.5;
    //                if (TmaxWork > 3000.0)
    //                    break;
    //                if (!TryEvalHint(TmaxWork, hint, out fmax))
    //                    return false;
    //            }
    //            else if (hint == Phases.Liquid)
    //            {
    //                TminWork = Math.Max(150.0, TminWork - 25.0);
    //                if (!TryEvalHint(TminWork, hint, out fmin))
    //                    return false;
    //            }
    //            else
    //            {
    //                TmaxWork *= 1.5;
    //                if (TmaxWork > 3000.0)
    //                    break;
    //                if (!TryEvalHint(TmaxWork, hint, out fmax))
    //                    return false;
    //            }
    //            expand++;
    //        }

    //        if (Math.Sign(fmin) == Math.Sign(fmax))
    //            return false;

    //        const int maxIts = 220;
    //        const double absTolH = 1e-3;
    //        const double relTolH = 1e-10;

    //        double a = TminWork, b = TmaxWork;
    //        double fa = fmin, fb = fmax;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double m = 0.5 * (a + b);
    //            if (!TryEvalHint(m, hint, out double fm))
    //                return false;

    //            if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
    //            {
    //                Tsol2 = m;
    //                return true;
    //            }

    //            if (Math.Sign(fm) == Math.Sign(fa))
    //            {
    //                a = m;
    //                fa = fm;
    //            }
    //            else
    //            {
    //                b = m;
    //                fb = fm;
    //            }
    //        }

    //        return false;
    //    }

    //    foreach (var hint in tries)
    //    {
    //        if (TrySolveHint(hint, out double Tsol3))
    //        {
    //            if (hint == Phases.Unknown)
    //                local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol3));
    //            else
    //                local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol3), hint);

    //            local.ClearTwoPhase();
    //            return;
    //        }
    //    }

    //    throw new InvalidOperationException(
    //        $"UpdatePH failed. Could not solve single-phase at P={pTarget.Pascal} Pa, h={hJ} J/kg. " +
    //        $"Tried hints: {string.Join(",", tries)} " +
    //        (sat != null ? $"(Tsat={TsatK}K, hL={hL}, hV={hV})" : "(no saturation info)"));
    //}

    public static void UpdatePH(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    {
        if (local == null)
            throw new ArgumentNullException(nameof(local));
        if (pTarget == null)
            throw new ArgumentNullException(nameof(pTarget));
        if (hTarget == null)
            throw new ArgumentNullException(nameof(hTarget));

        local.ClearTwoPhase();

        double Pc = local.Critical.Pressure.Pascal;
        double Tc = local.Critical.Temperature.Kelvin;
        double hJ = hTarget.JoulePerKilogram;

        // ---------------------------------------------------------------------
        // A) TWO-PHASE DETECTION (subcritical only)
        //
        // CoolProp-style structure: try to identify two-phase states first and set quality/state explicitly,
        // rather than accidentally solving them as single-phase. [1](https://cersonsky-lab.github.io/cbe710-notes/lecture_files/Lecture24.html)[2](https://github.com/CoolProp/CoolProp/discussions/2403)
        //
        // We try:
        //   1) EOS saturation solver SolveAtP(P)
        //   2) fallback: ancillary Tsat(P) + ancillary rhoL/rhoV(Tsat) -> construct SatResult record
        // ---------------------------------------------------------------------
        SaturationSolver.SatResult? sat = null;
        double TsatK = double.NaN;
        double hL = double.NaN, hV = double.NaN;

        if (pTarget.Pascal < Pc)
        {
            // 1) Try EOS saturation solver at P
            try
            {
                sat = local.SolveAtP(pTarget);
            }
            catch
            {
                sat = null;
            }

            // 2) Fallback saturation from ancillaries (no SolveAtT required)
            if (sat == null)
            {
                try
                {
                    var Tsat = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
                    var rhoLmol = LiquidDensity.CalculateDensity(Tsat);
                    var rhoVmol = VaporDensity.CalculateDensity(Tsat);

                    sat = new SaturationSolver.SatResult(Tsat, pTarget, rhoLmol, rhoVmol);
                }
                catch
                {
                    sat = null;
                }
            }

            // If saturation info exists, compute hL/hV from your EOS endpoints and check inside dome
            if (sat != null)
            {
                TsatK = sat.T.Kelvin;

                try
                {
                    // Saturated endpoints from your EOS at (Tsat, rho_sat)
                    var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * local.MolarMass };
                    var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * local.MolarMass };

                    hL = L.Enthalpy.SI;
                    hV = V.Enthalpy.SI;

                    double dh = hV - hL;

                    // Near critical dh collapses, quality ill-conditioned; treat as single-phase (CoolProp also treats this region specially) [3](https://colab.research.google.com/github/JMQuinlan/Thermo/blob/main/Thermo_Calc_Example.ipynb)
                    if (dh > 1e-3) // J/kg
                    {
                        // Buffer so CoolProp-based hMix doesn't "just miss" [hL,hV] due to small model differences
                        double epsH = Math.Max(50.0, 1e-6 * Math.Abs(dh));

                        if (hJ >= hL - epsH && hJ <= hV + epsH)
                        {
                            double q = (hJ - hL) / dh;
                            q = Math.Max(0.0, Math.Min(1.0, q));
                            local.SetTwoPhase(sat, q);
                            return;
                        }
                    }
                }
                catch
                {
                    // If endpoint evaluation fails, ignore saturation info and proceed single-phase
                    sat = null;
                }
            }
        }

        // ---------------------------------------------------------------------
        // B) SINGLE-PHASE FLASH (P,H -> T)
        //
        // 1) Blind solve first using your robust unhinted UpdatePT(P,T).
        //    This avoids dependency on saturation and is analogous to CoolProp flash routines using bounded solvers. [2](https://github.com/CoolProp/CoolProp/discussions/2403)[1](https://cersonsky-lab.github.io/cbe710-notes/lecture_files/Lecture24.html)
        //
        // 2) If blind fails, try hinted solves (Gas/Liquid/Unknown), with Tsat-side constraints if available.
        // ---------------------------------------------------------------------

        double Tmin = Math.Max(local.TripleLiquid.Temperature.Kelvin + 1e-3, 150.0);
        double Tmax = 2500.0;

        bool TryEvalBlind(double TK, out double f)
        {
            f = double.NaN;
            try
            {
                local.UpdatePT(pTarget, Temperature.FromKelvin(TK));
                f = local.Enthalpy.SI - hJ;
                return double.IsFinite(f);
            }
            catch
            {
                return false;
            }
        }

        bool TrySolveBlind(out double Tsol)
        {
            Tsol = double.NaN;

            double T0 = local.Temperature?.Kelvin ?? (Tc + 20.0);
            T0 = Math.Max(Tmin, Math.Min(Tmax, T0));

            // Find a valid start
            if (!TryEvalBlind(T0, out double f0))
            {
                double step = 20.0;
                bool found = false;

                for (int k = 0; k < 60; k++)
                {
                    double Tup = Math.Min(Tmax, T0 + step);
                    if (TryEvalBlind(Tup, out f0))
                    { T0 = Tup; found = true; break; }

                    double Tdn = Math.Max(Tmin, T0 - step);
                    if (TryEvalBlind(Tdn, out f0))
                    { T0 = Tdn; found = true; break; }

                    step *= 1.25;
                }

                if (!found)
                    return false;
            }

            // Bracket
            double Tlo = T0, Thi = T0;
            double flo = f0, fhi = f0;

            double step2 = 25.0;
            for (int k = 0; k < 140 && Math.Sign(flo) == Math.Sign(fhi); k++)
            {
                bool moved = false;

                double Tup = Math.Min(Tmax, Thi + step2);
                if (Tup > Thi + 1e-12 && TryEvalBlind(Tup, out double fup))
                {
                    Thi = Tup;
                    fhi = fup;
                    moved = true;
                }

                double Tdn = Math.Max(Tmin, Tlo - step2);
                if (Tdn < Tlo - 1e-12 && TryEvalBlind(Tdn, out double fdn))
                {
                    Tlo = Tdn;
                    flo = fdn;
                    moved = true;
                }

                if (!moved)
                    break;
                step2 *= 1.18;
            }

            if (Math.Sign(flo) == Math.Sign(fhi))
                return false;

            // Bisection
            const int maxIts = 240;
            const double absTolH = 1e-3;
            const double relTolH = 1e-10;

            for (int iter = 0; iter < maxIts; iter++)
            {
                double Tm = 0.5 * (Tlo + Thi);
                if (!TryEvalBlind(Tm, out double fm))
                    return false;

                if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
                {
                    Tsol = Tm;
                    return true;
                }

                if (Math.Sign(fm) == Math.Sign(flo))
                { Tlo = Tm; flo = fm; }
                else
                { Thi = Tm; fhi = fm; }
            }

            return false;
        }

        if (TrySolveBlind(out double Tblind))
        {
            local.UpdatePT(pTarget, Temperature.FromKelvin(Tblind));
            local.ClearTwoPhase();
            return;
        }

        // -------------------------
        // Hinted solve fallback
        // -------------------------
        bool TryEvalHint(double TK, Phases hint, out double f)
        {
            f = double.NaN;
            try
            {
                if (sat != null)
                {
                    if (hint == Phases.Liquid)
                        TK = Math.Min(TK, TsatK - 1e-6);
                    if (hint == Phases.Gas)
                        TK = Math.Max(TK, TsatK + 1e-6);
                }

                local.UpdatePT(pTarget, Temperature.FromKelvin(TK), hint);
                f = local.Enthalpy.SI - hJ;
                return double.IsFinite(f);
            }
            catch
            {
                return false;
            }
        }

        bool TrySolveHint(Phases hint, out double Tsol)
        {
            Tsol = double.NaN;

            double TminWork = Tmin;
            double TmaxWork = Tmax;

            if (sat != null)
            {
                if (hint == Phases.Liquid)
                    TmaxWork = Math.Min(TmaxWork, TsatK - 1e-6);
                if (hint == Phases.Gas)
                    TminWork = Math.Max(TminWork, TsatK + 1e-6);
            }

            if (TmaxWork <= TminWork + 1e-9)
                return false;

            if (!TryEvalHint(TminWork, hint, out double fmin))
                return false;
            if (!TryEvalHint(TmaxWork, hint, out double fmax))
                return false;

            int expand = 0;
            while (Math.Sign(fmin) == Math.Sign(fmax) && expand < 12)
            {
                if (hint == Phases.Gas)
                {
                    TmaxWork *= 1.5;
                    if (TmaxWork > 3000.0)
                        break;
                    if (!TryEvalHint(TmaxWork, hint, out fmax))
                        return false;
                }
                else if (hint == Phases.Liquid)
                {
                    TminWork = Math.Max(150.0, TminWork - 25.0);
                    if (!TryEvalHint(TminWork, hint, out fmin))
                        return false;
                }
                else
                {
                    TmaxWork *= 1.5;
                    if (TmaxWork > 3000.0)
                        break;
                    if (!TryEvalHint(TmaxWork, hint, out fmax))
                        return false;
                }
                expand++;
            }

            if (Math.Sign(fmin) == Math.Sign(fmax))
                return false;

            const int maxIts = 220;
            const double absTolH = 1e-3;
            const double relTolH = 1e-10;

            double a = TminWork, b = TmaxWork;
            double fa = fmin, fb = fmax;

            for (int iter = 0; iter < maxIts; iter++)
            {
                double m = 0.5 * (a + b);
                if (!TryEvalHint(m, hint, out double fm))
                    return false;

                if (Math.Abs(fm) < absTolH || Math.Abs(fm) / Math.Max(Math.Abs(hJ), 1.0) < relTolH)
                {
                    Tsol = m;
                    return true;
                }

                if (Math.Sign(fm) == Math.Sign(fa))
                { a = m; fa = fm; }
                else
                { b = m; fb = fm; }
            }

            return false;
        }

        // Decide hint order
        Phases[] hints = sat == null
            ? new[] { Phases.Unknown, Phases.Gas, Phases.Liquid }
            : (hJ < hL ? new[] { Phases.Liquid, Phases.Gas, Phases.Unknown }
                      : hJ > hV ? new[] { Phases.Gas, Phases.Liquid, Phases.Unknown }
                               : new[] { Phases.Unknown, Phases.Gas, Phases.Liquid });

        foreach (var hint in hints)
        {
            if (TrySolveHint(hint, out double Tsol))
            {
                if (hint == Phases.Unknown)
                    local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol));
                else
                    local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol), hint);

                local.ClearTwoPhase();
                return;
            }
        }

        throw new InvalidOperationException(
            $"UpdatePH failed. Could not solve single-phase at P={pTarget.Pascal} Pa, h={hJ} J/kg. " +
            $"Tried hints: {string.Join(",", hints)} " +
            (sat != null ? $"(Tsat={TsatK}K, hL={hL}, hV={hV})" : "(no saturation info)"));
    }

    public static void UpdatePX(this Ammonia local, Pressure pTarget, double quality)
    {
        if (local == null)
            throw new ArgumentNullException(nameof(local));
        if (pTarget == null)
            throw new ArgumentNullException(nameof(pTarget));
        if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
            throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

        if (pTarget.Pascal >= local.Critical.Pressure.Pascal)
            throw new InvalidOperationException($"UpdatePX invalid at/above Pc. P={pTarget.Pascal} Pa.");

        // ---- 1) Try EOS saturation solver ----
        SaturationSolver.SatResult? sat = null;
        try
        {
            sat = local.SolveAtP(pTarget);
        }
        catch
        {
            sat = null;
        }

        // ---- 2) Validate consistency of the sat result ----
        // We must ensure sat.T actually corresponds to pTarget.
        // Use your existing saturation correlation Psat(T) as a cheap consistency check.
        bool satOk = false;
        if (sat != null)
        {
            try
            {
                // correlation pressure at sat.T
                var psatFromT = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationPressure(sat.T);

                // Acceptable mismatch (tune if you want). We only need this as a sanity check.
                // If sat.T is wildly wrong, this check will fail by orders of magnitude.
                double absErr = Math.Abs(psatFromT.Pascal - pTarget.Pascal);
                satOk = absErr < 200.0; // 200 Pa sanity band
            }
            catch
            {
                satOk = false;
            }
        }

        // ---- 3) Fallback: Tsat(P) + rhoL/rhoV(Tsat) from ancillaries ----
        if (sat == null || !satOk)
        {
            var Tsat = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
            var rhoL = LiquidDensity.CalculateDensity(Tsat);
            var rhoV = VaporDensity.CalculateDensity(Tsat);
            sat = new SaturationSolver.SatResult(Tsat, pTarget, rhoL, rhoV);
        }
        else
        {
            // Ensure PQ contract: enforce exact pressure
            sat = sat with { Psat = pTarget };
        }

        local.SetTwoPhase(sat!, quality);
    }


    //public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t)
    //{
    //    if (local == null)
    //        throw new ArgumentNullException(nameof(local));
    //    if (pTarget == null)
    //        throw new ArgumentNullException(nameof(pTarget));
    //    if (t == null)
    //        throw new ArgumentNullException(nameof(t));

    //    local.ClearTwoPhase();
    //    local.Temperature = t;

    //    double Tc = local.Critical.Temperature.Kelvin;
    //    double Pc = local.Critical.Pressure.Pascal;
    //    double rhoRed = local.Critical.MolarDensity.MolesPerCubicMeter;

    //    bool subcritical = (pTarget.Pascal < Pc) && (t.Kelvin < Tc);

    //    // Ideal gas molar density guess
    //    double rhoIdeal = pTarget.Pascal / (local.GasConstant.SI * t.Kelvin);
    //    rhoIdeal = Math.Max(rhoIdeal, 1e-12);

    //    // Tsat(P): EOS if possible, else ancillary fallback
    //    Temperature TsatP = local.Critical.Temperature;
    //    SaturationSolver.SatResult? satP = null;

    //    if (subcritical)
    //    {
    //        try
    //        {
    //            satP = local.SolveAtP(pTarget);
    //            TsatP = satP.T;
    //        }
    //        catch
    //        {
    //            TsatP = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
    //            satP = null;
    //        }
    //    }

    //    bool wantLiquid = subcritical && (t.Kelvin < TsatP.Kelvin);

    //    // Ancillary density guess (only valid in [Ttriple, Tc])
    //    Molarity? rhoGuessMol = null;
    //    if (subcritical &&
    //        t.Kelvin > local.TripleLiquid.Temperature.Kelvin + 1e-6 &&
    //        t.Kelvin < Tc - 1e-6)
    //    {
    //        rhoGuessMol = wantLiquid ? LiquidDensity.CalculateDensity(t) : VaporDensity.CalculateDensity(t);
    //        if (rhoGuessMol != null && !double.IsFinite(rhoGuessMol.MolesPerCubicMeter))
    //            rhoGuessMol = null;
    //    }

    //    double PressureAt(double rhoMolLocal)
    //    {
    //        local.Density = Molarity.FromMolesPerCubicMeter(rhoMolLocal) * local.MolarMass;
    //        return local.Pressure.Pascal;
    //    }

    //    // Continuation Newton for vapor at high subcritical pressures:
    //    // keeps solution on low-density branch by limiting upward moves.
    //    bool TryVaporContinuation()
    //    {
    //        if (!subcritical)
    //            return false;
    //        if (wantLiquid)
    //            return false;
    //        if (!(pTarget.Pascal > 0.60 * Pc))
    //            return false; // only near-critical region

    //        // Use satP bounds if available to prevent crossing
    //        double rhoMax = 0.95 * rhoRed;
    //        if (satP != null)
    //        {
    //            double rhoL_satP = satP.RhomolarL.MolesPerCubicMeter;
    //            rhoMax = Math.Min(rhoMax, 0.80 * rhoL_satP);
    //        }

    //        // Start close to vapor: ideal or ancillary vapor guess
    //        double rho = rhoGuessMol?.MolesPerCubicMeter ?? rhoIdeal;
    //        rho = Math.Max(1e-12, Math.Min(rho, rhoMax));

    //        const int maxIts = 60;
    //        const double absTolPa = 0.5;
    //        const double relTol = 1e-10;

    //        for (int i = 0; i < maxIts; i++)
    //        {
    //            local.Density = Molarity.FromMolesPerCubicMeter(rho) * local.MolarMass;

    //            double f = local.Pressure.Pascal - pTarget.Pascal;

    //            if (Math.Abs(f) < absTolPa || Math.Abs(f) / Math.Max(pTarget.Pascal, 1.0) < relTol)
    //                return true;

    //            double dfdRho = local.dp_drhomolar_constT_SI;
    //            if (!(dfdRho > 0) || !double.IsFinite(dfdRho))
    //                return false;

    //            double step = -f / dfdRho;
    //            double rhoNew = rho + step;

    //            // Trust region: prevent jumping into dense-root basin
    //            // allow at most +50% increase per step, and at most -50% decrease per step
    //            double rhoUp = Math.Min(rho * 1.5, rhoMax);
    //            double rhoDn = Math.Max(rho / 1.5, 1e-12);
    //            rhoNew = Math.Max(rhoDn, Math.Min(rhoUp, rhoNew));

    //            // If we are stuck (no progress), expand trust region slightly
    //            if (Math.Abs(rhoNew - rho) / Math.Max(rho, 1.0) < 1e-12)
    //                rhoNew = Math.Min(rho * 1.05, rhoMax);

    //            rho = rhoNew;
    //        }

    //        return false;
    //    }

    //    bool TrySolve(bool liquidBranch)
    //    {
    //        // If vapor in near-critical region, try continuation first
    //        if (subcritical && !liquidBranch)
    //        {
    //            if (TryVaporContinuation())
    //                return true;
    //        }

    //        double rhoLoMol, rhoHiMol;
    //        double rhoHiCap = double.PositiveInfinity;

    //        if (subcritical)
    //        {
    //            if (liquidBranch)
    //            {
    //                double seed = rhoGuessMol?.MolesPerCubicMeter ?? Math.Max(2.0 * rhoRed, 50.0 * rhoIdeal);

    //                rhoLoMol = Math.Max(1e-12, 0.50 * seed);
    //                rhoHiMol = Math.Max(rhoLoMol * 1.01, 2.00 * seed);
    //                rhoHiMol = Math.Min(rhoHiMol, 20.0 * rhoRed);

    //                local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(seed, rhoLoMol, rhoHiMol)) * local.MolarMass;
    //            }
    //            else
    //            {
    //                // Vapor bracket (use satP if available)
    //                if (satP != null)
    //                {
    //                    double rhoV_satP = satP.RhomolarV.MolesPerCubicMeter;
    //                    double rhoL_satP = satP.RhomolarL.MolesPerCubicMeter;

    //                    double seed = rhoGuessMol?.MolesPerCubicMeter ?? Math.Min(rhoIdeal, rhoV_satP);

    //                    rhoLoMol = Math.Max(1e-12, Math.Min(seed, rhoV_satP) / 50.0);

    //                    rhoHiCap = Math.Min(5.0 * rhoV_satP, 0.50 * rhoL_satP);
    //                    rhoHiCap = Math.Max(rhoHiCap, rhoLoMol * 1.01);

    //                    rhoHiMol = Math.Min(Math.Max(2.0 * rhoV_satP, 20.0 * rhoIdeal), rhoHiCap);
    //                    rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);

    //                    local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(seed, rhoLoMol, rhoHiMol)) * local.MolarMass;
    //                }
    //                else
    //                {
    //                    double seed = rhoGuessMol?.MolesPerCubicMeter ?? rhoIdeal;

    //                    rhoLoMol = Math.Max(1e-12, seed / 1e6);
    //                    rhoHiMol = Math.Max(rhoLoMol * 1.01, 20.0 * rhoIdeal);

    //                    rhoHiCap = Math.Min(0.95 * rhoRed, 50.0 * rhoIdeal);
    //                    rhoHiMol = Math.Min(rhoHiMol, rhoHiCap);
    //                    rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);

    //                    local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(seed, rhoLoMol, rhoHiMol)) * local.MolarMass;
    //                }
    //            }
    //        }
    //        else
    //        {
    //            // Supercritical/single-root region
    //            double seed = rhoGuessMol?.MolesPerCubicMeter ?? rhoIdeal;

    //            rhoLoMol = Math.Max(1e-12, seed / 50.0);
    //            rhoHiMol = Math.Max(seed * 50.0, 5.0 * rhoRed);
    //            if (rhoHiMol <= rhoLoMol)
    //                rhoHiMol = rhoLoMol * 1.01;

    //            local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(seed, rhoLoMol, rhoHiMol)) * local.MolarMass;
    //        }

    //        double fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
    //        double fHi = PressureAt(rhoHiMol) - pTarget.Pascal;

    //        for (int k = 0; k < 80 && Math.Sign(fLo) == Math.Sign(fHi); k++)
    //        {
    //            rhoLoMol = Math.Max(1e-12, rhoLoMol * 0.5);

    //            if (subcritical && !liquidBranch)
    //            {
    //                rhoHiMol = Math.Min(rhoHiMol * 1.4, rhoHiCap);
    //                rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);

    //                fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
    //                fHi = PressureAt(rhoHiMol) - pTarget.Pascal;

    //                if (double.IsFinite(rhoHiCap) &&
    //                    Math.Abs(rhoHiMol - rhoHiCap) / Math.Max(rhoHiCap, 1.0) < 1e-12 &&
    //                    Math.Sign(fLo) == Math.Sign(fHi))
    //                    break;
    //            }
    //            else
    //            {
    //                rhoHiMol *= 2.0;
    //                fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
    //                fHi = PressureAt(rhoHiMol) - pTarget.Pascal;
    //            }
    //        }

    //        if (Math.Sign(fLo) == Math.Sign(fHi))
    //            return false;

    //        const int maxIts = 80;
    //        const double absTolPa = 0.5;
    //        const double relTol = 1e-10;

    //        double rhoMol = 0.5 * (rhoLoMol + rhoHiMol);
    //        local.Density = Molarity.FromMolesPerCubicMeter(rhoMol) * local.MolarMass;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double f = local.Pressure.Pascal - pTarget.Pascal;

    //            if (Math.Abs(f) < absTolPa || Math.Abs(f) / Math.Max(pTarget.Pascal, 1.0) < relTol)
    //                return true;

    //            double dfdRho = local.dp_drhomolar_constT_SI;
    //            bool newtonOK = dfdRho > 0 && double.IsFinite(dfdRho);

    //            double rhoNew = newtonOK ? (rhoMol - f / dfdRho) : double.NaN;
    //            if (!double.IsFinite(rhoNew) || rhoNew <= rhoLoMol || rhoNew >= rhoHiMol)
    //                rhoNew = 0.5 * (rhoLoMol + rhoHiMol);

    //            double fNew = PressureAt(rhoNew) - pTarget.Pascal;

    //            if (Math.Sign(fNew) == Math.Sign(fLo))
    //            {
    //                rhoLoMol = rhoNew;
    //                fLo = fNew;
    //            }
    //            else
    //            {
    //                rhoHiMol = rhoNew;
    //                fHi = fNew;
    //            }

    //            rhoMol = rhoNew;
    //        }

    //        return false;
    //    }

    //    if (TrySolve(wantLiquid))
    //        return;

    //    if (subcritical && TrySolve(!wantLiquid))
    //        return;

    //    throw new InvalidOperationException($"UpdatePT: could not bracket/solve. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
    //}



    //public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint)
    //{
    //    if (local == null)
    //        throw new ArgumentNullException(nameof(local));
    //    if (pTarget == null)
    //        throw new ArgumentNullException(nameof(pTarget));
    //    if (t == null)
    //        throw new ArgumentNullException(nameof(t));

    //    local.ClearTwoPhase();
    //    local.Temperature = t;

    //    double rhoRed = Ammonia.Critical.MolarDensity.MolesPerCubicMeter;

    //    // Ideal-gas molar density (always safe)
    //    double rhoIdeal = pTarget.Pascal / (Ammonia.GasConstant.SI * t.Kelvin);
    //    rhoIdeal = Math.Max(rhoIdeal, 1e-12);

    //    // Try ancillaries for density guess
    //    Molarity? rhoGuessMol = null;
    //    if (t.Kelvin > Ammonia.TripleLiquid.Temperature.Kelvin + 1e-6 &&
    //        t.Kelvin < Ammonia.Critical.Temperature.Kelvin - 1e-6)
    //    {
    //        rhoGuessMol = (phaseHint == Phases.Liquid)
    //            ? LiquidDensity.CalculateDensity(t)
    //            : VaporDensity.CalculateDensity(t);

    //        if (rhoGuessMol != null && !double.IsFinite(rhoGuessMol.MolesPerCubicMeter))
    //            rhoGuessMol = null;
    //    }

    //    double rhoSeed = rhoGuessMol?.MolesPerCubicMeter ?? rhoIdeal;

    //    double PressureAt(double rhoMolLocal)
    //    {
    //        local.Density = Molarity.FromMolesPerCubicMeter(rhoMolLocal) * Ammonia.MolarMass;
    //        return local.Pressure.Pascal;
    //    }

    //    bool TrySolveWithHint(Phases hint)
    //    {
    //        double rhoLoMol, rhoHiMol;

    //        if (hint == Phases.Gas)
    //        {
    //            // Vapor bracket: keep low to avoid dense root, but allow enough room
    //            rhoLoMol = Math.Max(1e-12, rhoSeed / 1e6);
    //            rhoHiMol = Math.Min(Math.Max(20.0 * rhoIdeal, 2.0 * rhoSeed), 0.95 * rhoRed);
    //            rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);

    //            local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(rhoSeed, rhoLoMol, rhoHiMol)) * Ammonia.MolarMass;
    //        }
    //        else if (hint == Phases.Liquid)
    //        {
    //            // Liquid bracket: high density
    //            rhoLoMol = Math.Max(1e-12, 0.5 * rhoSeed);
    //            rhoHiMol = Math.Max(rhoLoMol * 1.01, 2.0 * rhoSeed);
    //            rhoHiMol = Math.Min(rhoHiMol, 20.0 * rhoRed);

    //            local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(rhoSeed, rhoLoMol, rhoHiMol)) * Ammonia.MolarMass;
    //        }
    //        else
    //        {
    //            // Unknown/supercritical: wide bracket
    //            rhoLoMol = Math.Max(1e-12, rhoSeed / 50.0);
    //            rhoHiMol = Math.Max(rhoSeed * 50.0, 5.0 * rhoRed);
    //            rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);

    //            local.Density = Molarity.FromMolesPerCubicMeter(Math.Clamp(rhoSeed, rhoLoMol, rhoHiMol)) * Ammonia.MolarMass;
    //        }

    //        // Bracket check
    //        double fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
    //        double fHi = PressureAt(rhoHiMol) - pTarget.Pascal;

    //        for (int k = 0; k < 60 && Math.Sign(fLo) == Math.Sign(fHi); k++)
    //        {
    //            rhoLoMol = Math.Max(1e-12, rhoLoMol * 0.5);

    //            if (hint == Phases.Gas)
    //            {
    //                rhoHiMol = Math.Min(rhoHiMol * 1.4, 0.95 * rhoRed);
    //                rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);
    //            }
    //            else
    //            {
    //                rhoHiMol *= 2.0;
    //            }

    //            fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
    //            fHi = PressureAt(rhoHiMol) - pTarget.Pascal;
    //        }

    //        if (Math.Sign(fLo) == Math.Sign(fHi))
    //            return false;

    //        // Safeguarded Newton
    //        const int maxIts = 80;
    //        const double absTolPa = 0.5;
    //        const double relTol = 1e-10;

    //        double rhoMol = 0.5 * (rhoLoMol + rhoHiMol);
    //        local.Density = Molarity.FromMolesPerCubicMeter(rhoMol) * Ammonia.MolarMass;

    //        for (int iter = 0; iter < maxIts; iter++)
    //        {
    //            double pCalc = local.Pressure.Pascal;
    //            double f = pCalc - pTarget.Pascal;

    //            if (Math.Abs(f) < absTolPa || Math.Abs(f) / Math.Max(pTarget.Pascal, 1.0) < relTol)
    //                return true;

    //            double dfdRho = local.dp_drhomolar_constT_SI;
    //            bool newtonOK = dfdRho > 0 && double.IsFinite(dfdRho);

    //            double rhoNew = newtonOK ? (rhoMol - f / dfdRho) : double.NaN;
    //            if (!double.IsFinite(rhoNew) || rhoNew <= rhoLoMol || rhoNew >= rhoHiMol)
    //                rhoNew = 0.5 * (rhoLoMol + rhoHiMol);

    //            double fNew = PressureAt(rhoNew) - pTarget.Pascal;

    //            if (Math.Sign(fNew) == Math.Sign(fLo))
    //            {
    //                rhoLoMol = rhoNew;
    //                fLo = fNew;
    //            }
    //            else
    //            {
    //                rhoHiMol = rhoNew;
    //                fHi = fNew;
    //            }

    //            rhoMol = rhoNew;
    //        }

    //        return false;
    //    }

    //    // 1) Try requested hint first
    //    if (TrySolveWithHint(phaseHint))
    //        return;

    //    // 2) If it fails to bracket/solve, fall back to opposite phase (CoolProp-like “don’t die on bad phase hint”) [1](https://github.com/CoolProp/CoolProp/discussions/2403)[2](https://coolprop.org/_static/doxygen/html/class_cool_prop_1_1_equation_of_state.html)
    //    if (phaseHint == Phases.Gas && TrySolveWithHint(Phases.Liquid))
    //        return;

    //    if (phaseHint == Phases.Liquid && TrySolveWithHint(Phases.Gas))
    //        return;

    //    // 3) Last resort: wide bracket
    //    if (TrySolveWithHint(Phases.Unknown))
    //        return;

    //    throw new InvalidOperationException($"UpdatePT(phaseHint={phaseHint}) could not bracket/solve. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
    //}

    


    public static void UpdatePS(this Ammonia local, Pressure pTarget, SpecificEntropy sTarget)
    {
        if (local == null)
            throw new ArgumentNullException(nameof(local));
        if (pTarget == null)
            throw new ArgumentNullException(nameof(pTarget));
        if (sTarget == null)
            throw new ArgumentNullException(nameof(sTarget));

        local.ClearTwoPhase();

        double Pc = local.Critical.Pressure.Pascal;
        double Tc = local.Critical.Temperature.Kelvin;
        double sJ = sTarget.JoulePerKilogramKelvin;

        // ---------------------------------------------------------------------
        // A) Two-phase detection (subcritical only)
        // ---------------------------------------------------------------------
        SaturationSolver.SatResult? sat = null;
        double TsatK = double.NaN;
        double sL = double.NaN, sV = double.NaN;

        if (pTarget.Pascal < Pc)
        {
            // Try EOS saturation solver
            try
            {
                sat = local.SolveAtP(pTarget);
            }
            catch
            {
                sat = null;
            }

            // Fallback: ancillary Tsat(P) + rhoL/rhoV(Tsat)
            if (sat == null)
            {
                try
                {
                    var Tsat = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
                    var rhoLmol = LiquidDensity.CalculateDensity(Tsat);
                    var rhoVmol = VaporDensity.CalculateDensity(Tsat);
                    sat = new SaturationSolver.SatResult(Tsat, pTarget, rhoLmol, rhoVmol);
                }
                catch
                {
                    sat = null;
                }
            }

            if (sat != null)
            {
                TsatK = sat.T.Kelvin;

                try
                {
                    var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * local.MolarMass };
                    var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * local.MolarMass };

                    sL = L.Entropy.SI; // J/kg/K
                    sV = V.Entropy.SI;

                    double ds = sV - sL;

                    // Near critical ds collapses; quality becomes ill-conditioned -> treat as single-phase
                    if (ds > 1e-6)
                    {
                        // Slop band to avoid falling through due to tiny endpoint mismatches
                        double epsS = Math.Max(1e-3, 1e-6 * Math.Abs(ds)); // J/kg/K

                        if (sJ >= sL - epsS && sJ <= sV + epsS)
                        {
                            double q = (sJ - sL) / ds;
                            q = Math.Max(0.0, Math.Min(1.0, q));
                            local.SetTwoPhase(sat, q);
                            return;
                        }
                    }
                }
                catch
                {
                    sat = null; // fall through to single-phase
                }
            }
        }

        // ---------------------------------------------------------------------
        // B) Single-phase solve for T: s(P,T) = sTarget
        // ---------------------------------------------------------------------
        double Tmin = Math.Max(local.TripleLiquid.Temperature.Kelvin + 1e-3, 150.0);
        double Tmax = 2500.0;

        // Pick preferred phase hint if sat endpoints exist
        Phases preferredHint = Phases.Unknown;
        if (sat != null)
        {
            if (sJ < sL)
                preferredHint = Phases.Liquid;
            else if (sJ > sV)
                preferredHint = Phases.Gas;
            else
                preferredHint = Phases.Unknown;
        }
        else
        {
            preferredHint = Phases.Unknown;
        }

        Phases[] hints = preferredHint == Phases.Liquid
            ? new[] { Phases.Liquid, Phases.Gas, Phases.Unknown }
            : preferredHint == Phases.Gas
                ? new[] { Phases.Gas, Phases.Liquid, Phases.Unknown }
                : new[] { Phases.Unknown, Phases.Gas, Phases.Liquid };

        bool TryEvalBlind(double TK, out double f)
        {
            f = double.NaN;
            try
            {
                local.UpdatePT(pTarget, Temperature.FromKelvin(TK));
                f = local.Entropy.SI - sJ;
                return double.IsFinite(f);
            }
            catch { return false; }
        }

        bool TryEvalHint(double TK, Phases hint, out double f)
        {
            f = double.NaN;
            try
            {
                if (sat != null)
                {
                    if (hint == Phases.Liquid)
                        TK = Math.Min(TK, TsatK - 1e-6);
                    if (hint == Phases.Gas)
                        TK = Math.Max(TK, TsatK + 1e-6);
                }

                if (hint == Phases.Unknown)
                    local.UpdatePT(pTarget, Temperature.FromKelvin(TK));
                else
                    local.UpdatePT(pTarget, Temperature.FromKelvin(TK), hint);

                f = local.Entropy.SI - sJ;
                return double.IsFinite(f);
            }
            catch { return false; }
        }

        bool SolveWith(Func<double, (bool ok, double f)> eval, out double Tsol)
        {
            Tsol = double.NaN;

            // Find a valid start near Tc
            double T0 = local.Temperature?.Kelvin ?? (Tc + 20.0);
            T0 = Math.Max(Tmin, Math.Min(Tmax, T0));

            bool ok0;
            double f0;
            (ok0, f0) = eval(T0);

            if (!ok0)
            {
                double step = 20.0;
                bool found = false;
                for (int k = 0; k < 60; k++)
                {
                    double Tup = Math.Min(Tmax, T0 + step);
                    (ok0, f0) = eval(Tup);
                    if (ok0)
                    { T0 = Tup; found = true; break; }

                    double Tdn = Math.Max(Tmin, T0 - step);
                    (ok0, f0) = eval(Tdn);
                    if (ok0)
                    { T0 = Tdn; found = true; break; }

                    step *= 1.25;
                }
                if (!found)
                    return false;
            }

            // Bracket
            double Tlo = T0, Thi = T0;
            double flo = f0, fhi = f0;

            double step2 = 25.0;
            for (int k = 0; k < 140 && Math.Sign(flo) == Math.Sign(fhi); k++)
            {
                bool moved = false;

                double Tup = Math.Min(Tmax, Thi + step2);
                var (oku, fu) = eval(Tup);
                if (oku)
                { Thi = Tup; fhi = fu; moved = true; }

                double Tdn = Math.Max(Tmin, Tlo - step2);
                var (okd, fd) = eval(Tdn);
                if (okd)
                { Tlo = Tdn; flo = fd; moved = true; }

                if (!moved)
                    break;
                step2 *= 1.18;
            }

            if (Math.Sign(flo) == Math.Sign(fhi))
                return false;

            // Bisection
            const int maxIts = 240;
            const double absTolS = 1e-4;   // J/kg/K
            const double relTolS = 1e-10;

            for (int iter = 0; iter < maxIts; iter++)
            {
                double Tm = 0.5 * (Tlo + Thi);
                var (okm, fm) = eval(Tm);
                if (!okm)
                    return false;

                if (Math.Abs(fm) < absTolS || Math.Abs(fm) / Math.Max(Math.Abs(sJ), 1.0) < relTolS)
                {
                    Tsol = Tm;
                    return true;
                }

                if (Math.Sign(fm) == Math.Sign(flo))
                { Tlo = Tm; flo = fm; }
                else
                { Thi = Tm; fhi = fm; }
            }

            return false;
        }

        // 1) Try blind solve first (no saturation dependency)
        if (SolveWith(TK => (TryEvalBlind(TK, out var f), f), out double Tblind))
        {
            local.UpdatePT(pTarget, Temperature.FromKelvin(Tblind));
            local.ClearTwoPhase();
            return;
        }

        // 2) Try hinted solves
        foreach (var hint in hints)
        {
            if (SolveWith(TK => (TryEvalHint(TK, hint, out var f), f), out double Tsol))
            {
                if (hint == Phases.Unknown)
                    local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol));
                else
                    local.UpdatePT(pTarget, Temperature.FromKelvin(Tsol), hint);

                local.ClearTwoPhase();
                return;
            }
        }

        throw new InvalidOperationException($"UpdatePS failed at P={pTarget.Pascal} Pa, s={sJ} J/kg/K.");
    }


        /// <summary>
        /// Update state from Temperature and vapor quality (mass quality in [0,1]).
        /// This is a two-phase flash for pure fluids below critical temperature.
        /// </summary>
        public static void UpdateTX(this EngineeringFluids.Fluids.Ammonia local, Temperature t, double quality)
        {
            if (local == null)
                throw new ArgumentNullException(nameof(local));
            if (t == null)
                throw new ArgumentNullException(nameof(t));
            if (double.IsNaN(quality) || quality < 0.0 || quality > 1.0)
                throw new ArgumentOutOfRangeException(nameof(quality), "Quality must be in [0,1].");

            // No two-phase exists at/above critical temperature
            if (t.Kelvin >= local.Critical.Temperature.Kelvin)
                throw new InvalidOperationException($"UpdateTX invalid at/above Tc. T={t.Kelvin} K.");

            // Guard very low temps (below triple)
            if (t.Kelvin < local.TripleLiquid.Temperature.Kelvin)
                throw new InvalidOperationException($"UpdateTX invalid below triple temperature. T={t.Kelvin} K.");

            // Clear and set two-phase
            local.ClearTwoPhase();

            // Saturation pressure at T (ancillary/correlation)
            Pressure psat = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationPressure(t);

            // Saturated densities at T (ancillaries)
            Molarity rhoL = LiquidDensity.CalculateDensity(t);
            Molarity rhoV = VaporDensity.CalculateDensity(t);

            // Build SatResult cache (your record type)
            var sat = new EngineeringFluids.Helmholtz.Solvers.SaturationSolver.SatResult(
                T: t,
                Psat: psat,
                RhomolarL: rhoL,
                RhomolarV: rhoV
            );

            local.SetTwoPhase(sat, quality);
        }

        /// <summary>
        /// Alias to match SharpFluids naming (Quality, Temperature).
        /// </summary>
        public static void UpdateXT(this EngineeringFluids.Fluids.Ammonia local, double quality, Temperature t)
            => UpdateTX(local, t, quality);
}

