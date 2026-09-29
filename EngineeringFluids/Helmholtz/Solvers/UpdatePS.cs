using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Pressure + entropy flash. Same two-part shape as UpdatePH.cs -
    // see its class remarks for the full rationale (two-phase ancillary check, warm-started
    // outer temperature root-find via Illinois with a safe-bisection fallback for the bracket-
    // collapse failure mode already found and fixed there). Entropy is monotonic increasing in
    // T at fixed P away from phase instabilities (dS/dT|P = Cp/T > 0), same as enthalpy
    // (dH/dT|P = Cp > 0), so the identical bracketing/root-find shape applies unchanged with
    // Entropy in place of Enthalpy at each trial T.
    public static void UpdatePS(this Ammonia local, Pressure pTarget, SpecificEntropy sTarget)
    {
        if (!double.IsFinite(pTarget.Pascal) || pTarget <= Pressure.Zero)
            throw new ArgumentOutOfRangeException(nameof(pTarget), $"UpdatePS: pressure must be positive. P={pTarget.Pascal} Pa.");
        if (!sTarget.HasValue())
            throw new ArgumentOutOfRangeException(nameof(sTarget), "UpdatePS: entropy must be finite.");

        local.ClearTwoPhase();

        Pressure Pc = local.Critical.Pressure;
        Temperature Tc = local.Critical.Temperature;
        Temperature Ttriple = local.TripleLiquid.Temperature;
        MolarMass M = local.MolarMass;

        Temperature Tmin = Ttriple + FlashTripleMargin;
        Temperature Tmax = FlashTmax;

        bool subcritical = pTarget < Pc;
        Temperature Tsat = Temperature.FromKelvin(double.NaN);
        SpecificEntropy sL = SpecificEntropy.FromJoulePerKilogramKelvin(double.NaN), sV = SpecificEntropy.FromJoulePerKilogramKelvin(double.NaN);
        Molarity rhoLGuess, rhoVGuess;
        bool sawDome = false;

        if (subcritical)
        {
            Tsat = SaturationTemperature.Temperature(pTarget);

            if (Tc - Tsat >= CriticalMargin)
            {
                rhoLGuess = BubbleDensity.Density(Tsat);
                rhoVGuess = DewDensity.Density(Tsat);

                // sL/sV come from the full EOS at Tsat/rho, the same way SetTwoPhase's own
                // SatLiquidState/SatVaporState will evaluate them afterward - not a separate
                // ancillary entropy fit. See UpdatePH.cs's hL/hV remarks for why this
                // self-consistency matters near the endpoints of the quality inversion.
                sL = new Ammonia { Temperature = Tsat, Density = rhoLGuess * M }.Entropy;
                sV = new Ammonia { Temperature = Tsat, Density = rhoVGuess * M }.Entropy;
                SpecificEntropy ds = sV - sL;

                if (ds > SpecificEntropy.FromJoulePerKilogramKelvin(1e-6))
                {
                    sawDome = true;

                    SpecificEntropy epsS = SpecificEntropy.Max(SpecificEntropy.FromJoulePerKilogramKelvin(1e-3), 1e-6 * ds.Abs());
                    if (sTarget >= sL - epsS && sTarget <= sV + epsS)
                    {
                        double q = Math.Clamp((double)((sTarget - sL) / ds), 0.0, 1.0);
                        var sat = new SaturationSolver.SatResult(Tsat, pTarget, rhoLGuess, rhoVGuess);
                        local.SetTwoPhase(sat, q);
                        return;
                    }
                }
            }
        }

        // ------------------------------------------------------------------
        // Single-phase flash: find T such that s(T) == sTarget at fixed P (see class remarks).
        // ------------------------------------------------------------------

        Phases fixedPhase = Phases.Unknown;
        Temperature lo, hi;

        if (sawDome && sTarget < sL)
        {
            fixedPhase = Phases.Liquid;
            lo = Tmin;
            hi = Tsat - FlashDomeMargin;
        }
        else if (sawDome) // sTarget > sV (the [sL,sV]+margin case already returned above)
        {
            fixedPhase = Phases.Gas;
            lo = Tsat + FlashDomeMargin;
            hi = Tmax;
        }
        else
        {
            lo = Tmin;
            hi = Tmax;
        }

        Molarity rhoGuessMolar = Molarity.Zero;
        bool haveGuess = false;

        SpecificEntropy EvalWarm(Temperature TK)
        {
            Phases p = fixedPhase;
            if (p == Phases.Unknown)
            {
                p = DeterminePhaseCheapPure(local, TK, pTarget);
                if (p == Phases.Twophase)
                {
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

            rhoGuessMolar = rhomolar;
            haveGuess = true;
            local.Temperature = TK;
            local.Density = rhomolar * M;
            return local.Entropy;
        }

        SpecificEntropy flo = EvalWarm(lo) - sTarget;
        SpecificEntropy fhi = EvalWarm(hi) - sTarget;

        if (Math.Sign(flo.SI) == Math.Sign(fhi.SI))
            throw new InvalidOperationException(
                $"UpdatePS: could not bracket a single-phase solution in [{lo.Kelvin:G6} K, {hi.Kelvin:G6} K]. " +
                $"P={pTarget.Pascal} Pa, s={sTarget.JoulePerKilogramKelvin} J/kg/K is likely outside the supported range.");

        // Illinois (regula falsi variant) - see UpdatePH.cs for why convergence is judged
        // only by the actual residual, never by bracket width, and why a safe-bisection
        // fallback is needed for the collapse failure mode.
        const int maxIts = 80;
        const double tolRel = 1e-9;
        SpecificEntropy sScale = SpecificEntropy.Max(sTarget.Abs(), SpecificEntropy.FromJoulePerKilogramKelvin(1.0));
        int side = 0;
        Temperature a = lo, b = hi;
        SpecificEntropy fa = flo, fb = fhi;
        Temperature lastC = Temperature.FromKelvin(double.NaN);
        SpecificEntropy lastFc = SpecificEntropy.FromJoulePerKilogramKelvin(double.NaN);

        for (int i = 0; i < maxIts; i++)
        {
            Temperature c = (fa * b - fb * a) / (fa - fb);
            c = c.Clamp(Temperature.Min(a, b), Temperature.Max(a, b));
            SpecificEntropy fc = EvalWarm(c) - sTarget;
            lastC = c;
            lastFc = fc;

            if (fc.Abs() <= tolRel * sScale)
                return;

            if (Math.Sign(fc.SI) == Math.Sign(fb.SI))
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

            if ((b - a).Abs() < FlashBracketCollapse)
            {
                SolveBySafeBisection(lo, hi, flo, fhi);
                return;
            }
        }

        if (lastFc.Abs() > 1e-6 * sScale)
        {
            SolveBySafeBisection(lo, hi, flo, fhi);
            return;
        }

        void SolveBySafeBisection(Temperature bA, Temperature bB, SpecificEntropy fA, SpecificEntropy fB)
        {
            for (int i = 0; i < 200; i++)
            {
                Temperature m = 0.5 * (bA + bB);
                SpecificEntropy fm = EvalWarm(m) - sTarget;

                if (fm.Abs() <= tolRel * sScale || (bB - bA).Abs() < FlashBisectionWidth)
                    return;

                if (Math.Sign(fm.SI) == Math.Sign(fA.SI))
                { bA = m; fA = fm; }
                else
                { bB = m; fB = fm; }
            }

            throw new InvalidOperationException(
                $"UpdatePS: single-phase root-find did not converge. P={pTarget.Pascal} Pa, s={sTarget.JoulePerKilogramKelvin} J/kg/K.");
        }
    }
}
