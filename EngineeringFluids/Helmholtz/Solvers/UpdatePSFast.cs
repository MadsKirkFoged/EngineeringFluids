using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using System;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{
    // Fast (AmmoniaDouble) pressure + entropy flash. Same two-part shape as UpdatePHFast.cs -
    // see its class remarks for the full rationale (two-phase ancillary check, warm-started
    // outer temperature root-find via Illinois with a safe-bisection fallback for the bracket-
    // collapse failure mode already found and fixed there). Entropy is monotonic increasing in
    // T at fixed P away from phase instabilities (dS/dT|P = Cp/T > 0), same as enthalpy
    // (dH/dT|P = Cp > 0), so the identical bracketing/root-find shape applies unchanged with
    // Entropy in place of Enthalpy at each trial T.
    public static void UpdatePS(this AmmoniaDouble local, double pTarget, double sTarget)
    {
        if (!double.IsFinite(pTarget) || pTarget <= 0)
            throw new ArgumentOutOfRangeException(nameof(pTarget), $"UpdatePS: pressure must be positive. P={pTarget} Pa.");
        if (!double.IsFinite(sTarget))
            throw new ArgumentOutOfRangeException(nameof(sTarget), "UpdatePS: entropy must be finite.");

        local.ClearTwoPhase();

        double Pc = local.Critical.Pressure;
        double Tc = local.Critical.Temperature;
        double Ttriple = local.TripleLiquid.Temperature;
        double M = local.MolarMass;

        double Tmin = Ttriple + 1e-3;
        const double Tmax = 2500.0;

        bool subcritical = pTarget < Pc;
        double Tsat = double.NaN, sL = double.NaN, sV = double.NaN;
        double rhoLGuess = double.NaN, rhoVGuess = double.NaN;
        bool sawDome = false;

        if (subcritical)
        {
            Tsat = SaturationTemperatureFast.Temperature((float)pTarget);

            if (Tc - Tsat >= 1.0)
            {
                rhoLGuess = BubbleDensityFast.Density((float)Tsat);
                rhoVGuess = DewDensityFast.Density((float)Tsat);

                // sL/sV come from the full EOS at Tsat/rho, the same way SetTwoPhase's own
                // SatLiquidState/SatVaporState will evaluate them afterward - not a separate
                // ancillary entropy fit. See UpdatePHFast.cs's hL/hV remarks for why this
                // self-consistency matters near the endpoints of the quality inversion.
                sL = new AmmoniaDouble { Temperature = Tsat, Density = rhoLGuess * M }.Entropy;
                sV = new AmmoniaDouble { Temperature = Tsat, Density = rhoVGuess * M }.Entropy;
                double ds = sV - sL;

                if (ds > 1e-6)
                {
                    sawDome = true;

                    double epsS = Math.Max(1e-3, 1e-6 * Math.Abs(ds));
                    if (sTarget >= sL - epsS && sTarget <= sV + epsS)
                    {
                        double q = Math.Clamp((sTarget - sL) / ds, 0.0, 1.0);
                        var sat = new SaturationSolver.SatResultDouble(Tsat, pTarget, rhoLGuess, rhoVGuess);
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
        double lo, hi;

        if (sawDome && sTarget < sL)
        {
            fixedPhase = Phases.Liquid;
            lo = Tmin;
            hi = Tsat - 1e-4;
        }
        else if (sawDome) // sTarget > sV (the [sL,sV]+margin case already returned above)
        {
            fixedPhase = Phases.Gas;
            lo = Tsat + 1e-4;
            hi = Tmax;
        }
        else
        {
            lo = Tmin;
            hi = Tmax;
        }

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

            rhoGuessMolar = rhomolar;
            haveGuess = true;
            local.Temperature = TK;
            local.Density = rhomolar * M;
            return local.Entropy;
        }

        double flo = EvalWarm(lo) - sTarget;
        double fhi = EvalWarm(hi) - sTarget;

        if (Math.Sign(flo) == Math.Sign(fhi))
            throw new InvalidOperationException(
                $"UpdatePS: could not bracket a single-phase solution in [{lo:G6} K, {hi:G6} K]. " +
                $"P={pTarget} Pa, s={sTarget} J/kg/K is likely outside the supported range.");

        // Illinois (regula falsi variant) - see UpdatePHFast.cs for why convergence is judged
        // only by the actual residual, never by bracket width, and why a safe-bisection
        // fallback is needed for the collapse failure mode.
        const int maxIts = 80;
        const double tolRel = 1e-9;
        int side = 0;
        double a = lo, b = hi, fa = flo, fb = fhi;
        double lastC = double.NaN, lastFc = double.NaN;

        for (int i = 0; i < maxIts; i++)
        {
            double c = (fa * b - fb * a) / (fa - fb);
            c = Math.Clamp(c, Math.Min(a, b), Math.Max(a, b));
            double fc = EvalWarm(c) - sTarget;
            lastC = c;
            lastFc = fc;

            if (Math.Abs(fc) <= tolRel * Math.Max(Math.Abs(sTarget), 1.0))
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

            if (Math.Abs(b - a) < 1e-9)
            {
                SolveBySafeBisection(lo, hi, flo, fhi);
                return;
            }
        }

        if (Math.Abs(lastFc) > 1e-6 * Math.Max(Math.Abs(sTarget), 1.0))
        {
            SolveBySafeBisection(lo, hi, flo, fhi);
            return;
        }

        void SolveBySafeBisection(double bA, double bB, double fA, double fB)
        {
            for (int i = 0; i < 200; i++)
            {
                double m = 0.5 * (bA + bB);
                double fm = EvalWarm(m) - sTarget;

                if (Math.Abs(fm) <= tolRel * Math.Max(Math.Abs(sTarget), 1.0) || Math.Abs(bB - bA) < 1e-12)
                    return;

                if (Math.Sign(fm) == Math.Sign(fA))
                { bA = m; fA = fm; }
                else
                { bB = m; fB = fm; }
            }

            throw new InvalidOperationException(
                $"UpdatePS: single-phase root-find did not converge. P={pTarget} Pa, s={sTarget} J/kg/K.");
        }
    }
}
