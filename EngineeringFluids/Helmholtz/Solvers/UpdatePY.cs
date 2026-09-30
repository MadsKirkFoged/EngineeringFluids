using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Temperature search limits shared by the single-phase (P,H) and (P,S) flashes.
    private static readonly Temperature FlashTmax = Temperature.FromKelvin(2500.0);
    private static readonly Temperature FlashTripleMargin = Temperature.FromKelvin(1e-3);
    private static readonly Temperature FlashDomeMargin = Temperature.FromKelvin(1e-4);
    private static readonly Temperature FlashTwoPhaseNudge = Temperature.FromKelvin(1e-3);
    private static readonly Temperature FlashBracketCollapse = Temperature.FromKelvin(1e-9);
    private static readonly Temperature FlashBisectionWidth = Temperature.FromKelvin(1e-12);

    // The exact (CoolProp-compatible) flashes search up to 1.5 x the EOS's maximum temperature, as CoolProp's
    // own P+H/P+S flash does (HSU_P_flash), so both reject the same out-of-range inputs.
    private static readonly Temperature ExactFlashTmax = Temperature.FromKelvin(1.5 * 725.0);

    // Once a bisection has pinned the root's temperature to 1e-12 K, a residual this small is the inner density
    // solve's noise (near the critical point h and s are extremely sensitive to density), not a wrong state.
    private const double FlashNoiseTolRel = 1e-6;

    // The exact flashes treat a target within this fraction of the dome width from a saturation line as ON that
    // line (quality 0 or 1). Exact saturation still differs from CoolProp's by ~1e-10, so a state CoolProp puts
    // exactly on the saturation line would otherwise land a hair outside the dome - or its bracket would collapse.
    private const double DomeEdgeBand = 1e-9;

    private enum FlashProperty { Enthalpy, Entropy }

    // Pressure + enthalpy flash. Unlike UpdatePT/UpdatePX, there is no
    // direct ancillary shortcut for (P,H): the result can be single-phase OR two-phase
    // depending on where hTarget falls, so this has to combine the cheap two-phase check from
    // UpdatePX.cs with an actual root-find over temperature for the single-phase case.
    //
    // A) Two-phase check: same ancillaries as UpdatePX (Tsat/rhoL/rhoV/hL/hV), reused here to
    //    see whether hTarget falls inside the saturation dome at Tsat(P). If so, quality is a
    //    direct linear interpolation - no iteration needed, same O(1) cost as UpdatePX.
    // B) Single-phase flash: false position (Illinois variant) over T at fixed P, bounded to
    //    the correct side of the dome (when known). Illinois converges superlinearly, so this
    //    typically needs only a handful of iterations, versus 200+ for plain bisection. Each
    //    iteration calls the internal
    //    density solver (SolveRhoMolar_TP, from UpdatePT.cs - private but visible here
    //    since this is the same partial class) directly with a density guess WARM-STARTED
    //    from the previous trial temperature, rather than going through the public UpdatePT
    //    (which re-detects phase and rebuilds a guess from ancillaries from scratch on every
    //    call). Density varies smoothly with T at fixed P, so once the outer search is
    //    anywhere near the root, consecutive trial temperatures are close together and the
    //    previous solution is an excellent starting point - this collapses the inner density
    //    solve to 1-2 Newton steps per outer iteration instead of rebuilding it from an
    //    ancillary guess every time. The net win still scales with how many outer iterations
    //    the bracket needs: a handful for a tight bracket near the dome, more for a bracket
    //    spanning a wide T range (e.g. a far-superheated gas state bracketed by [Tsat, 2500K]).
    //    [benchmark-guided]
    //
    // In the last millikelvin below Tc the two-phase band collapses to near-zero width and even
    // the exact densities are ill-conditioned (see UpdatePX) - this falls through to the
    // single-phase search instead of trying to resolve a vanishing dome.
    //
    // Separately, at LOW absolute pressure with a two-phase result very close to q=0 or q=1,
    // the reported density can be less precise in relative terms than elsewhere: rhoV is tiny
    // there, so d(specific volume)/dq ~ 1/rhoV is huge, and any residual mismatch between this
    // EOS's hL/hV and CoolProp's true values (however small) gets amplified through the
    // quality-from-enthalpy inversion. This is inherent to inverting h back to a quality at
    // low P/low q - not present in UpdatePX, which is given quality directly. The absolute
    // density error stays tiny even so; only points below ~1 MPa with quality very close to an
    // endpoint are affected. [benchmark-guided]
    //
    // A target within 50 J/kg of either end of the dome counts as two-phase (see the buffer below), so a
    // state just outside the dome comes back as saturated liquid/vapor. UpdatePHExact has no such buffer.
    public static void UpdatePH(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    {
        ValidatePY(pTarget, hTarget.JoulePerKilogram, FlashProperty.Enthalpy, nameof(hTarget), nameof(UpdatePH));
        UpdatePY(local, pTarget, hTarget.JoulePerKilogram, FlashProperty.Enthalpy, exact: false, nameof(UpdatePH));
    }

    // Pressure + entropy flash. Same two-part shape as UpdatePH - see its remarks for the full rationale
    // (two-phase ancillary check, warm-started outer temperature root-find via Illinois with a safe-bisection
    // fallback for the bracket-collapse failure mode already found and fixed there). Entropy is monotonic
    // increasing in T at fixed P away from phase instabilities (dS/dT|P = Cp/T > 0), same as enthalpy
    // (dH/dT|P = Cp > 0), so the identical bracketing/root-find shape applies unchanged.
    public static void UpdatePS(this Ammonia local, Pressure pTarget, SpecificEntropy sTarget)
    {
        ValidatePY(pTarget, sTarget.JoulePerKilogramKelvin, FlashProperty.Entropy, nameof(sTarget), nameof(UpdatePS));
        UpdatePY(local, pTarget, sTarget.JoulePerKilogramKelvin, FlashProperty.Entropy, exact: false, nameof(UpdatePS));
    }

    // Like UpdatePH, but the two-phase part is solved against the EOS (SolveAtP) instead of read from the
    // ancillary fits, with no buffer around the dome and CoolProp's temperature limits - so it matches CoolProp
    // (and SharpFluids) to ~1e-10 in two-phase too, right up to the critical point. Costs a few SolveAtT calls
    // more than UpdatePH when the pressure is below Pc.
    public static void UpdatePHExact(this Ammonia local, Pressure pTarget, Enthalpy hTarget)
    {
        ValidatePY(pTarget, hTarget.JoulePerKilogram, FlashProperty.Enthalpy, nameof(hTarget), nameof(UpdatePHExact));
        UpdatePY(local, pTarget, hTarget.JoulePerKilogram, FlashProperty.Enthalpy, exact: true, nameof(UpdatePHExact));
    }

    // Like UpdatePS, with the two-phase part solved exactly - see UpdatePHExact.
    public static void UpdatePSExact(this Ammonia local, Pressure pTarget, SpecificEntropy sTarget)
    {
        ValidatePY(pTarget, sTarget.JoulePerKilogramKelvin, FlashProperty.Entropy, nameof(sTarget), nameof(UpdatePSExact));
        UpdatePY(local, pTarget, sTarget.JoulePerKilogramKelvin, FlashProperty.Entropy, exact: true, nameof(UpdatePSExact));
    }

    private static void ValidatePY(Pressure pTarget, double yTarget, FlashProperty y, string yName, string name)
    {
        if (!double.IsFinite(pTarget.Pascal) || pTarget <= Pressure.Zero)
            throw new ArgumentOutOfRangeException(nameof(pTarget), $"{name}: pressure must be positive. P={pTarget.Pascal} Pa.");
        if (!double.IsFinite(yTarget))
            throw new ArgumentOutOfRangeException(yName, $"{name}: {(y == FlashProperty.Enthalpy ? "enthalpy" : "entropy")} must be finite.");
    }

    // The property being matched, in SI (J/kg or J/kg/K). Plain doubles because one solver serves both.
    private static double Read(Ammonia a, FlashProperty y)
        => y == FlashProperty.Enthalpy ? a.Enthalpy.JoulePerKilogram : a.Entropy.JoulePerKilogramKelvin;

    private static string Describe(FlashProperty y, double value)
        => y == FlashProperty.Enthalpy ? $"h={value} J/kg" : $"s={value} J/kg/K";

    private static void UpdatePY(Ammonia local, Pressure pTarget, double yTarget, FlashProperty y, bool exact, string name)
    {
        local.ClearTwoPhase();

        Pressure Pc = local.Critical.Pressure;
        Temperature Tc = local.Critical.Temperature;
        Temperature Ttriple = local.TripleLiquid.Temperature;
        MolarMass M = local.MolarMass;

        Temperature Tmin = exact ? Ttriple : Ttriple + FlashTripleMargin;
        Temperature Tmax = exact ? ExactFlashTmax : FlashTmax;
        Temperature domeMargin = exact ? Temperature.Zero : FlashDomeMargin;

        // Narrowest dome still treated as one, and the buffer around it (fast mode only)
        double minDomeWidth = y == FlashProperty.Enthalpy ? 1e-3 : 1e-6;
        double bufferAbs = y == FlashProperty.Enthalpy ? 50.0 : 1e-3;

        bool subcritical = pTarget < Pc;
        Temperature Tsat = Temperature.FromKelvin(double.NaN);
        double yL = double.NaN, yV = double.NaN;
        Ammonia? satL = null, satV = null;
        bool sawDome = false;

        if (subcritical)
        {
            SaturationSolver.SatResult? sat = exact ? TrySolveAtP(local, pTarget) : AncillarySaturation(local, pTarget);

            if (sat is not null)
            {
                Tsat = sat.T;

                // yL/yV must come from the SAME source SetTwoPhase's own SatLiquidState/
                // SatVaporState will use afterward (the full EOS evaluated at Tsat/rhoL/rhoV),
                // not separate bubble/dew enthalpy polynomial fits. Those
                // fits carry their own small, independent approximation error relative to the
                // EOS - fine for UpdatePX (which only needs a self-consistent quality/rho pair
                // once), but here that mismatch was showing up amplified into density error
                // near q~0: the reciprocal specific-volume mixing rule is very sensitive to a
                // small quality error there (1/rhoV >> 1/rhoL). Evaluating yL/yV the same way
                // the final state will be reported keeps this self-consistent. [benchmark-guided]
                satL = new Ammonia { Temperature = Tsat, Density = sat.RhomolarL * M };
                satV = new Ammonia { Temperature = Tsat, Density = sat.RhomolarV * M };
                yL = Read(satL, y);
                yV = Read(satV, y);
                double dy = yV - yL;

                if (dy > minDomeWidth)
                {
                    sawDome = true;

                    // Fast mode: small buffer so a target right at the endpoint isn't missed due to tiny
                    // model differences between the ancillary fit and the EOS itself.
                    double buffer = exact ? DomeEdgeBand * dy : Math.Max(bufferAbs, 1e-6 * Math.Abs(dy));
                    if (yTarget >= yL - buffer && yTarget <= yV + buffer)
                    {
                        double q = Math.Clamp((yTarget - yL) / dy, 0.0, 1.0);
                        // the saturated states just evaluated are the ones the two-phase properties read
                        local.SetTwoPhase(new SaturationSolver.SatResult(Tsat, pTarget, sat.RhomolarL, sat.RhomolarV), q, Phases.Twophase, satL, satV);
                        return;
                    }
                }
            }
        }

        // ------------------------------------------------------------------
        // Single-phase flash: find T such that y(T) == yTarget at fixed P (see UpdatePH remarks).
        // ------------------------------------------------------------------

        Phases fixedPhase = Phases.Unknown;
        Temperature lo, hi;

        if (sawDome && yTarget < yL)
        {
            fixedPhase = Phases.Liquid;
            lo = Tmin;
            hi = Tsat - domeMargin;
        }
        else if (sawDome) // yTarget > yV (the [yL,yV]+buffer case already returned above)
        {
            fixedPhase = Phases.Gas;
            lo = Tsat + domeMargin;
            hi = Tmax;
        }
        else
        {
            // Either supercritical (P >= Pc) or the near-critical dome-skip case: phase can
            // legitimately change across the search range (Liquid for T<Tc, Supercritical for
            // T>=Tc, per this EOS's own phase convention), so it's re-checked per trial inside
            // EvalWarm rather than fixed once up front.
            lo = Tmin;
            hi = Tmax;
        }

        // Warm-started inner evaluator: calls the density solver directly with the PREVIOUS
        // trial's converged density as the starting guess, instead of going through the public
        // UpdatePT (which always rebuilds a guess from ancillaries) - EXCEPT when the new trial
        // temperature is far from the previous one, where a stale warm-start guess is actively
        // dangerous rather than merely suboptimal. [benchmark-guided: the very first two trials
        // are the two ends of the WHOLE bracket - e.g. Tmin near the triple point vs. Tsat
        // possibly 100+ K higher for the liquid branch - the worst possible pair to reuse a
        // guess between (liquid density there can differ by several-fold). Warm-starting only
        // pays off once Illinois has narrowed the bracket enough that consecutive trials are
        // genuinely close - but "genuinely close" turned out to have no safe fixed threshold:
        // liquid density is sensitive enough to temperature that even ~1% T changes have been
        // seen to carry a warm-start guess outside where Newton can recover (rhoMin is a firm
        // 1.05x critical density there). Rather than chase a threshold, this tries the warm
        // start and falls back to the always-reliable ancillary guess if it fails - the same
        // safeguarded-fallback shape as the density Newton solve itself.]
        //
        // A warm start can also CONVERGE to the wrong root: started from a low supercritical density, the liquid
        // solve found a density inside the EOS's unstable loop where the pressure matches too (248 instead of
        // 474 kg/m3 at 121.55 bar, 375.5 K). That silently wrong h misled the outer search, so warm-started
        // liquid roots are checked against the saturated-liquid density and redone cold when implausible.
        Molarity rhoGuessMolar = Molarity.Zero;
        bool haveGuess = false;

        double EvalWarm(Temperature TK)
        {
            Phases p = fixedPhase;
            if (p == Phases.Unknown)
            {
                p = DeterminePhaseCheapPure(local, TK, pTarget);
                if (p == Phases.Twophase)
                {
                    // Only reachable in the near-critical-skip case, landing in the razor-thin
                    // two-phase ambiguity band UpdatePT itself would also reject - nudge away.
                    TK += FlashTwoPhaseNudge;
                    p = DeterminePhaseCheapPure(local, TK, pTarget);
                }
            }

            Molarity rhomolar;
            if (haveGuess)
            {
                try
                {
                    rhomolar = SolveRhoMolar_TP(local, TK, pTarget, p, rhoGuessMolar);
                    if (!IsPlausibleRoot(local, TK, p, rhomolar))
                        rhomolar = SolveRhoMolar_TP(local, TK, pTarget, p, GuessRhoMolar(local, TK, pTarget, p, out _));
                }
                catch (InvalidOperationException)
                {
                    rhomolar = SolveRhoMolar_TP(local, TK, pTarget, p, GuessRhoMolar(local, TK, pTarget, p, out _));
                }
            }
            else
            {
                rhomolar = SolveRhoMolar_TP(local, TK, pTarget, p, GuessRhoMolar(local, TK, pTarget, p, out _));
            }

            rhoGuessMolar = rhomolar; // warm-start candidate for the next trial
            haveGuess = true;
            local.Temperature = TK;
            local.Density = rhomolar * M;
            return Read(local, y);
        }

        // Newton first, when the dome fixes the phase: dy/dT at fixed P is Cp for enthalpy and Cp/T for entropy,
        // both from the EOS evaluation the property read just cached, and the saturated state's Cp gives a start a
        // few kelvin from the root. The dome side of the bracket is known without evaluating it (y is monotonic in
        // T and y(Tsat) is yL/yV), so steps are kept inside [lo, hi] and bisect when they would leave it. Any
        // failure to converge falls through to the Illinois search below, which is unchanged.
        if (sawDome && satL is not null && satV is not null)
        {
            try
            {
                if (TryNewtonT(fixedPhase == Phases.Gas ? satV : satL, fixedPhase == Phases.Gas ? yV : yL))
                    return;
            }
            catch (InvalidOperationException)
            {
            }
            haveGuess = false;
        }

        bool TryNewtonT(Ammonia saturated, double ySat)
        {
            const double tolNewton = 1e-9; // same residual test as the Illinois search
            double scale = Math.Max(Math.Abs(yTarget), 1.0);
            double cpSat = saturated.Cp.JoulePerKilogramKelvin;
            double dySat = yTarget - ySat;
            Temperature T = y == FlashProperty.Enthalpy
                ? Tsat + Temperature.FromKelvin(dySat / cpSat)
                : Tsat * Math.Exp(dySat / cpSat);

            Temperature a = lo, b = hi; // y(a) < yTarget < y(b)
            if (!(T > a && T < b))
                T = 0.5 * (a + b);

            for (int i = 0; i < 25; i++)
            {
                double f = EvalWarm(T) - yTarget;
                if (Math.Abs(f) <= tolNewton * scale)
                    return true;

                if (f < 0) a = T; else b = T;
                if ((b - a).Abs() < FlashBracketCollapse)
                    return false;

                double cp = local.Cp.JoulePerKilogramKelvin;
                double dydT = y == FlashProperty.Enthalpy ? cp : cp / T.Kelvin;
                Temperature next = T - Temperature.FromKelvin(f / dydT);
                T = next > a && next < b ? next : 0.5 * (a + b);
            }
            return false;
        }

        double flo = EvalWarm(lo) - yTarget;
        double fhi = EvalWarm(hi) - yTarget;

        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"{name}: could not bracket a single-phase solution in [{lo.Kelvin:G6} K, {hi.Kelvin:G6} K]. " +
                $"P={pTarget.Pascal} Pa, {Describe(y, yTarget)} is likely outside the supported range.");

        // Illinois (regula falsi variant): robust - the bracket is never lost - and converges
        // superlinearly for a smooth, monotonic function like enthalpy vs. temperature at
        // fixed pressure, unlike plain bisection or unguarded secant/Newton.
        //
        // Convergence is judged ONLY by the actual residual, never by bracket width: the
        // Illinois modification below (halving the stale side's f-value to stop the bracket
        // stalling) deliberately drives that side's TRACKED f toward zero, which pulls the
        // next c toward it - so the bracket can collapse to near-zero width from that
        // artificial decay alone, at a point whose TRUE residual is still far from zero.
        // [benchmark-guided: this "|b-a| small => converged" shortcut was accepting a state
        // with h off by ~5% at some high-pressure liquid points - the pressure at that state
        // matched the target essentially exactly, since the density solve at each trial T is
        // independently correct; only the OUTER root-find had stopped too early.]
        const int maxIts = 80;
        const double tolRel = 1e-9;
        double yScale = Math.Max(Math.Abs(yTarget), 1.0);
        int side = 0;
        Temperature a = lo, b = hi;
        double fa = flo, fb = fhi;
        double lastFc = double.NaN;

        for (int i = 0; i < maxIts; i++)
        {
            Temperature c = (fa * b - fb * a) / (fa - fb);
            c = c.Clamp(Temperature.Min(a, b), Temperature.Max(a, b));
            double fc = EvalWarm(c) - yTarget; // side effect: local is now at T=c
            lastFc = fc;

            if (Math.Abs(fc) <= tolRel * yScale)
                return;

            if (Math.Sign(fc) == Math.Sign(fb))
            {
                b = c; fb = fc;
                if (side == -1) fa *= 0.5;
                side = -1;
            }
            else
            {
                a = c; fa = fc;
                if (side == 1) fb *= 0.5;
                side = 1;
            }

            // The bracket has collapsed (from the Illinois modification's artificial decay,
            // per the remarks above) without the residual actually converging - a genuine
            // stall, not success. Restart the search as a plain bisection over the ORIGINAL
            // bracket, which has no artificial f-value decay to cause this failure mode
            // (slower per step, but its "small interval => converged" reasoning is actually
            // valid, and it always finds the true root of a correctly-bracketed monotonic
            // function).
            if ((b - a).Abs() < FlashBracketCollapse)
            {
                SolveBySafeBisection(lo, hi, flo, fhi);
                return;
            }
        }

        // Illinois normally converges in well under maxIts for a smooth monotonic function;
        // if it didn't, fall back to safe bisection rather than accept an unconverged result.
        if (Math.Abs(lastFc) > 1e-6 * yScale)
        {
            SolveBySafeBisection(lo, hi, flo, fhi);
            return;
        }

        void SolveBySafeBisection(Temperature bA, Temperature bB, double fA, double fB)
        {
            for (int i = 0; i < 200; i++)
            {
                Temperature m = 0.5 * (bA + bB);
                double fm = EvalWarm(m) - yTarget;

                if (Math.Abs(fm) <= tolRel * yScale || (bB - bA).Abs() < FlashBisectionWidth)
                {
                    // A pinned bracket with a large residual straddles a jump in y(T), not a root:
                    // report it instead of returning a state that does not match the input.
                    if (!(Math.Abs(fm) <= FlashNoiseTolRel * yScale))
                        break;
                    return;
                }

                if (Math.Sign(fm) == Math.Sign(fA))
                { bA = m; fA = fm; }
                else
                { bB = m; fB = fm; }
            }

            throw new InvalidOperationException(
                $"{name}: single-phase root-find did not converge. P={pTarget.Pascal} Pa, {Describe(y, yTarget)}.");
        }
    }

    // Fast mode's saturation: the ancillary tables, and none within the last millikelvin below Tc (see UpdatePX)
    private static SaturationSolver.SatResult? AncillarySaturation(Ammonia local, Pressure p)
    {
        Temperature Tsat = SaturationTemperature.Temperature(p);
        if (local.Critical.Temperature - Tsat < CriticalMargin)
            return null;
        return new SaturationSolver.SatResult(Tsat, p, BubbleDensity.Density(Tsat), DewDensity.Density(Tsat));
    }

    // Exact mode's saturation: solved against the EOS, which converges to within 6 microK of Tc (tested); if it fails
    // closer than that, the dome is too narrow to matter and the single-phase search covers it, as in fast mode.
    private static SaturationSolver.SatResult? TrySolveAtP(Ammonia local, Pressure p)
    {
        try
        {
            return local.SolveAtP(p);
        }
        catch (InvalidOperationException)
        {
            return null;
        }
    }

    // A stable liquid at p > psat(T) is always at least as dense as the saturated liquid, so a liquid-branch
    // root well below the bubble density is a spurious root of the EOS's unstable loop.
    private static bool IsPlausibleRoot(Ammonia a, Temperature T, Phases phase, Molarity rhomolar)
    {
        if (phase != Phases.Liquid || T < a.TripleLiquid.Temperature || T >= a.Critical.Temperature - TemperatureMargin)
            return true;
        return rhomolar >= (1 - 1e-3) * BubbleDensity.Density(T);
    }
}
