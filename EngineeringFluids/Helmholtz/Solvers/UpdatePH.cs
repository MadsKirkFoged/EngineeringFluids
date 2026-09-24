using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Fast (AmmoniaDouble) pressure + enthalpy flash. Unlike UpdatePT/UpdatePX, there is no
    // direct ancillary shortcut for (P,H): the result can be single-phase OR two-phase
    // depending on where hTarget falls, so this has to combine the cheap two-phase check from
    // UpdatePXFast.cs with an actual root-find over temperature for the single-phase case.
    //
    // A) Two-phase check: same ancillaries as UpdatePX (Tsat/rhoL/rhoV/hL/hV), reused here to
    //    see whether hTarget falls inside the saturation dome at Tsat(P). If so, quality is a
    //    direct linear interpolation - no iteration needed, same O(1) cost as UpdatePX.
    // B) Single-phase flash: false position (Illinois variant) over T at fixed P, bounded to
    //    the correct side of the dome (when known). Illinois converges superlinearly, so this
    //    typically needs only a handful of iterations, versus 200+ for plain bisection. Each
    //    iteration calls the internal
    //    density solver (SolveRhoMolar_TP, from UpdatePTFast.cs - private but visible here
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
    // Near the critical point the two-phase band collapses to near-zero width, and both the
    // ancillaries and CoolProp's own saturation solver lose reliability resolving it (same
    // fragility already found for UpdatePT/UpdatePX) - this falls through to the single-phase
    // search instead of trying to resolve a vanishing dome. [benchmark-guided]
    //
    // Separately, at LOW absolute pressure with a two-phase result very close to q=0 or q=1,
    // the reported density can be less precise in relative terms than elsewhere: rhoV is tiny
    // there, so d(specific volume)/dq ~ 1/rhoV is huge, and any residual mismatch between this
    // EOS's hL/hV and CoolProp's true values (however small) gets amplified through the
    // quality-from-enthalpy inversion. This is inherent to inverting h back to a quality at
    // low P/low q - not present in UpdatePX, which is given quality directly. The absolute
    // density error stays tiny even so; only points below ~1 MPa with quality very close to an
    // endpoint are affected. [benchmark-guided]
    public static void UpdatePH(this AmmoniaDouble local, double pTarget, double hTarget)
    {
        if (!double.IsFinite(pTarget) || pTarget <= 0)
            throw new ArgumentOutOfRangeException(nameof(pTarget), $"UpdatePH: pressure must be positive. P={pTarget} Pa.");
        if (!double.IsFinite(hTarget))
            throw new ArgumentOutOfRangeException(nameof(hTarget), "UpdatePH: enthalpy must be finite.");

        local.ClearTwoPhase();

        double Pc = local.Critical.Pressure;
        double Tc = local.Critical.Temperature;
        double Ttriple = local.TripleLiquid.Temperature;
        double M = local.MolarMass;
        double rhoRed = local.Critical.MolarDensity;

        double Tmin = Ttriple + 1e-3;
        const double Tmax = 2500.0;

        bool subcritical = pTarget < Pc;
        double Tsat = double.NaN, hL = double.NaN, hV = double.NaN;
        double rhoLGuess = double.NaN, rhoVGuess = double.NaN;
        bool sawDome = false;

        if (subcritical)
        {
            Tsat = SaturationTemperatureFast.Temperature((float)pTarget);

            if (Tc - Tsat >= 1.0)
            {
                rhoLGuess = BubbleDensityFast.Density((float)Tsat);
                rhoVGuess = DewDensityFast.Density((float)Tsat);

                // hL/hV must come from the SAME source SetTwoPhase's own SatLiquidState/
                // SatVaporState will use afterward (the full EOS evaluated at Tsat/rhoL/rhoV),
                // not separate bubble/dew enthalpy polynomial fits. Those
                // fits carry their own small, independent approximation error relative to the
                // EOS - fine for UpdatePX (which only needs a self-consistent quality/rho pair
                // once), but here that mismatch was showing up amplified into density error
                // near q~0: the reciprocal specific-volume mixing rule is very sensitive to a
                // small quality error there (1/rhoV >> 1/rhoL). Evaluating hL/hV the same way
                // the final state will be reported keeps this self-consistent. [benchmark-guided]
                hL = new AmmoniaDouble { Temperature = Tsat, Density = rhoLGuess * M }.Enthalpy;
                hV = new AmmoniaDouble { Temperature = Tsat, Density = rhoVGuess * M }.Enthalpy;
                double dh = hV - hL;

                if (dh > 1e-3)
                {
                    sawDome = true;

                    // Small buffer so a target right at the endpoint isn't missed due to tiny
                    // model differences between the ancillary fit and the EOS itself.
                    double epsH = Math.Max(50.0, 1e-6 * Math.Abs(dh));
                    if (hTarget >= hL - epsH && hTarget <= hV + epsH)
                    {
                        double q = Math.Clamp((hTarget - hL) / dh, 0.0, 1.0);
                        var sat = new SaturationSolver.SatResultDouble(Tsat, pTarget, rhoLGuess, rhoVGuess);
                        local.SetTwoPhase(sat, q);
                        return;
                    }
                }
            }
        }

        // ------------------------------------------------------------------
        // Single-phase flash: find T such that h(T) == hTarget at fixed P (see class remarks).
        // ------------------------------------------------------------------

        Phases fixedPhase = Phases.Unknown;
        double lo, hi;

        if (sawDome && hTarget < hL)
        {
            fixedPhase = Phases.Liquid;
            lo = Tmin;
            hi = Tsat - 1e-4;
        }
        else if (sawDome) // hTarget > hV (the [hL,hV]+margin case already returned above)
        {
            fixedPhase = Phases.Gas;
            lo = Tsat + 1e-4;
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
        double rhoGuessMolar = 0.0;
        bool haveGuess = false;

        double EvalWarm(double TK)
        {
            Phases p = fixedPhase;
            if (p == Phases.Unknown)
            {
                p = DeterminePhaseCheapPure(local, TK, pTarget);
                if (p == Phases.Twophase)
                {
                    // Only reachable in the near-critical-skip case, landing in the razor-thin
                    // two-phase ambiguity band UpdatePT itself would also reject - nudge away.
                    TK += 1e-3;
                    p = DeterminePhaseCheapPure(local, TK, pTarget);
                }
            }

            double rhomolar;
            if (haveGuess)
            {
                try
                {
                    rhomolar = SolveRhoMolar_TP(local, TK, pTarget, p, rhoGuessMolar);
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
            return local.Enthalpy;
        }

        double flo = EvalWarm(lo) - hTarget;
        double fhi = EvalWarm(hi) - hTarget;

        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"UpdatePH: could not bracket a single-phase solution in [{lo:G6} K, {hi:G6} K]. " +
                $"P={pTarget} Pa, h={hTarget} J/kg is likely outside the supported range.");

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
        int side = 0;
        double a = lo, b = hi, fa = flo, fb = fhi;
        double lastC = double.NaN, lastFc = double.NaN;

        for (int i = 0; i < maxIts; i++)
        {
            double c = (fa * b - fb * a) / (fa - fb);
            c = Math.Clamp(c, Math.Min(a, b), Math.Max(a, b));
            double fc = EvalWarm(c) - hTarget; // side effect: local is now at T=c
            lastC = c;
            lastFc = fc;

            if (Math.Abs(fc) <= tolRel * Math.Max(Math.Abs(hTarget), 1.0))
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
            if (Math.Abs(b - a) < 1e-9)
            {
                SolveBySafeBisection(lo, hi, flo, fhi);
                return;
            }
        }

        // Illinois normally converges in well under maxIts for a smooth monotonic function;
        // if it didn't, fall back to safe bisection rather than accept an unconverged result.
        if (Math.Abs(lastFc) > 1e-6 * Math.Max(Math.Abs(hTarget), 1.0))
        {
            SolveBySafeBisection(lo, hi, flo, fhi);
            return;
        }

        void SolveBySafeBisection(double bA, double bB, double fA, double fB)
        {
            for (int i = 0; i < 200; i++)
            {
                double m = 0.5 * (bA + bB);
                double fm = EvalWarm(m) - hTarget;

                if (Math.Abs(fm) <= tolRel * Math.Max(Math.Abs(hTarget), 1.0) || Math.Abs(bB - bA) < 1e-12)
                    return;

                if (Math.Sign(fm) == Math.Sign(fA))
                { bA = m; fA = fm; }
                else
                { bB = m; fB = fm; }
            }

            throw new InvalidOperationException(
                $"UpdatePH: single-phase root-find did not converge. P={pTarget} Pa, h={hTarget} J/kg.");
        }
    }
}
