using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Ancillary;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits;
using static EngineeringFluids.Helmholtz.Phase;

public static partial class Update
{




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

//    public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint)
//    {
//        if (local == null)
//            throw new ArgumentNullException(nameof(local));
//        if (pTarget == null)
//            throw new ArgumentNullException(nameof(pTarget));
//        if (t == null)
//            throw new ArgumentNullException(nameof(t));

//        local.ClearTwoPhase();
//        local.Temperature = t;

//        double rhoRed = local.Critical.MolarDensity.MolesPerCubicMeter;

//        // Ideal gas molar density (safe default)
//        double rhoIdeal = pTarget.Pascal / (local.GasConstant.SI * t.Kelvin);
//        rhoIdeal = Math.Max(rhoIdeal, 1e-12);

//        // Ancillary density guess if available
//        Molarity? rhoGuessMol = null;
//        if (t.Kelvin > local.TripleLiquid.Temperature.Kelvin + 1e-6 &&
//            t.Kelvin < local.Critical.Temperature.Kelvin - 1e-6)
//        {
//            rhoGuessMol = (phaseHint == Phases.Liquid)
//                ? LiquidDensity.CalculateDensity(t)
//                : VaporDensity.CalculateDensity(t);

//            if (rhoGuessMol != null && !double.IsFinite(rhoGuessMol.MolesPerCubicMeter))
//                rhoGuessMol = null;
//        }

//        double rhoSeed = rhoGuessMol?.MolesPerCubicMeter ?? rhoIdeal;

//        double PressureAt(double rhoMolLocal)
//        {
//            local.Density = Molarity.FromMolesPerCubicMeter(rhoMolLocal) * local.MolarMass;
//            return local.Pressure.Pascal;
//        }

//        // Build bracket by phase
//        double rhoLoMol, rhoHiMol;

//        if (phaseHint == Phases.Gas)
//        {
//            rhoLoMol = Math.Max(1e-12, rhoSeed / 1e6);
//            rhoHiMol = Math.Min(Math.Max(20.0 * rhoIdeal, 2.0 * rhoSeed), 0.95 * rhoRed);
//            rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);
//        }
//        else if (phaseHint == Phases.Liquid)
//        {
//            rhoLoMol = Math.Max(1e-12, 0.5 * rhoSeed);
//            rhoHiMol = Math.Max(rhoLoMol * 1.01, 2.0 * rhoSeed);
//            rhoHiMol = Math.Min(rhoHiMol, 20.0 * rhoRed);
//        }
//        else
//        {
//            rhoLoMol = Math.Max(1e-12, rhoSeed / 50.0);
//            rhoHiMol = Math.Max(rhoSeed * 50.0, 5.0 * rhoRed);
//            rhoHiMol = Math.Max(rhoHiMol, rhoLoMol * 1.01);
//        }

//        // Bracket check
//        double fLo = PressureAt(rhoLoMol) - pTarget.Pascal;
//        double fHi = PressureAt(rhoHiMol) - pTarget.Pascal;

//        for (int k = 0; k < 60 && Math.Sign(fLo) == Math.Sign(fHi); k++)
//        {
//            rhoLoMol = Math.Max(1e-12, rhoLoMol * 0.5);

//            if (phaseHint == Phases.Gas)
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
//            throw new InvalidOperationException($"UpdatePT(phaseHint={phaseHint}) could not bracket. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");

//        // Safeguarded Newton
//        const int maxIts = 80;
//        const double absTolPa = 0.5;
//        const double relTol = 1e-10;

//        double rhoMol = 0.5 * (rhoLoMol + rhoHiMol);
//        local.Density = Molarity.FromMolesPerCubicMeter(rhoMol) * local.MolarMass;

//        for (int iter = 0; iter < maxIts; iter++)
//        {
//            double pCalc = local.Pressure.Pascal;
//            double f = pCalc - pTarget.Pascal;

//            if (Math.Abs(f) < absTolPa || Math.Abs(f) / Math.Max(pTarget.Pascal, 1.0) < relTol)
//                return;

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

//        throw new InvalidOperationException($"UpdatePT(phaseHint={phaseHint}) did not converge. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
//    }



        // Public API stays the same
        //public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t)
        //    => UpdatePTCore(local, pTarget, t, Phases.Unknown, strictHint: false);

        //public static void UpdatePT(this Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint)
        //    => UpdatePTCore(local, pTarget, t, phaseHint, strictHint: true);


    //private static void UpdatePTCore(Ammonia local, Pressure pTarget, Temperature t, Phases phaseHint, bool strictHint)
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

    //    // Ideal gas molar density guess (always safe)
    //    double rhoIdeal = pTarget.Pascal / (local.GasConstant.SI * t.Kelvin);
    //    rhoIdeal = Math.Max(rhoIdeal, 1e-12);

    //    // Tsat(P): used only in AUTO mode to pick branch + provide vapor caps
    //    Temperature TsatP = local.Critical.Temperature;
    //    SaturationSolver.SatResult? satP = null;

    //    bool wantLiquid;

    //    if (strictHint)
    //    {
    //        if (phaseHint != Phases.Liquid && phaseHint != Phases.Gas)
    //            throw new ArgumentException($"UpdatePT strict requires phaseHint=Liquid or Gas, got {phaseHint}");

    //        wantLiquid = (phaseHint == Phases.Liquid);
    //    }
    //    else
    //    {
    //        if (subcritical)
    //        {
    //            try
    //            {
    //                satP = local.SolveAtP(pTarget);
    //                TsatP = satP.T;
    //            }
    //            catch
    //            {
    //                TsatP = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationTemperature(pTarget);
    //                satP = null;
    //            }
    //        }

    //        wantLiquid = subcritical && (t.Kelvin < TsatP.Kelvin);
    //    }

    //    // Ancillary guess only valid in (Ttriple, Tc)
    //    Molarity? rhoGuessMol = null;
    //    if (subcritical &&
    //        t.Kelvin > local.TripleLiquid.Temperature.Kelvin + 1e-6 &&
    //        t.Kelvin < Tc - 1e-6)
    //    {
    //        bool guessLiquid = strictHint ? (phaseHint == Phases.Liquid) : wantLiquid;
    //        rhoGuessMol = guessLiquid ? LiquidDensity.CalculateDensity(t) : VaporDensity.CalculateDensity(t);

    //        if (rhoGuessMol != null && !double.IsFinite(rhoGuessMol.MolesPerCubicMeter))
    //            rhoGuessMol = null;
    //    }

    //    double PressureAt(double rhoMolLocal)
    //    {
    //        local.Density = Molarity.FromMolesPerCubicMeter(rhoMolLocal) * local.MolarMass;
    //        return local.Pressure.Pascal;
    //    }

    //    // ------------------------------------------------------------
    //    // Near-critical vapor continuation (only in AUTO vapor solves)
    //    // Keeps solution on low-density branch via trust region.
    //    // ------------------------------------------------------------
    //    bool TryVaporContinuation()
    //    {
    //        if (strictHint)
    //            return false;
    //        if (!subcritical)
    //            return false;
    //        if (wantLiquid)
    //            return false;
    //        if (!(pTarget.Pascal > 0.60 * Pc))
    //            return false;

    //        // ensure satP if possible (better caps)
    //        if (satP == null)
    //        {
    //            try
    //            { satP = local.SolveAtP(pTarget); }
    //            catch { satP = null; }
    //        }

    //        double rhoMax = 0.95 * rhoRed;
    //        if (satP != null)
    //        {
    //            double rhoL_satP = satP.RhomolarL.MolesPerCubicMeter;
    //            rhoMax = Math.Min(rhoMax, 0.80 * rhoL_satP);
    //        }

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

    //            double rhoNew = rho + (-f / dfdRho);

    //            // Trust region
    //            double rhoUp = Math.Min(rho * 1.5, rhoMax);
    //            double rhoDn = Math.Max(rho / 1.5, 1e-12);
    //            rhoNew = Math.Max(rhoDn, Math.Min(rhoUp, rhoNew));

    //            if (Math.Abs(rhoNew - rho) / Math.Max(rho, 1.0) < 1e-12)
    //                rhoNew = Math.Min(rho * 1.05, rhoMax);

    //            rho = rhoNew;
    //        }

    //        return false;
    //    }

    //    // ------------------------------------------------------------
    //    // Solve one branch (liquidBranch=true/false)
    //    // ------------------------------------------------------------
    //    bool TrySolveBranch(bool liquidBranch)
    //    {
    //        if (!strictHint && subcritical && !liquidBranch)
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
    //                // Vapor bracket with sat bounds if available (AUTO may have satP; strict vapor probably won't)
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

    //        // Safeguarded Newton
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

    //    // ------------------------------------------------------------
    //    // Solve order
    //    // ------------------------------------------------------------
    //    if (strictHint)
    //    {
    //        bool isLiquid = (phaseHint == Phases.Liquid);
    //        if (TrySolveBranch(isLiquid))
    //            return;

    //        throw new InvalidOperationException(
    //            $"UpdatePT(phaseHint={phaseHint}) could not bracket/solve. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
    //    }
    //    else
    //    {
    //        if (TrySolveBranch(wantLiquid))
    //            return;

    //        if (subcritical && TrySolveBranch(!wantLiquid))
    //            return;

    //        throw new InvalidOperationException(
    //            $"UpdatePT: could not bracket/solve. T={t.Kelvin} K, P={pTarget.Pascal} Pa.");
    //    }
    //}


    




}