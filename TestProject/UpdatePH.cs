using EngineeringFluids.Fluids;
using EngineeringUnits;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;

namespace TestProject;

[TestClass]
public class CoolPropOracle_UpdatePH_FastTests
{
    private const double Tc = 405.56;
    private const double Pc = 11363400.0;
    private const double Ptriple = 6091.223108650368;

    private static double RelErr(double expected, double actual)
    {
        double denom = Math.Max(Math.Abs(expected), 1.0);
        return Math.Abs(actual - expected) / denom;
    }

    // ------------------------------------------------------------------
    // Fixed-grid single-phase: build a known (P,T) state via SharpFluids/CoolProp, feed its
    // enthalpy back into our UpdatePH, and check it recovers the same T/rho/s/u.
    // ------------------------------------------------------------------
    public static IEnumerable<object[]> SinglePhasePoints()
    {
        // Gas, far from the dome
        yield return new object[] { 2_100_000.0, 373.15 };
        // Gas, near the dome
        yield return new object[] { 5_630_000.0, 373.15 };
        // Liquid
        yield return new object[] { 2_100_000.0, 280.0 };
        yield return new object[] { 8_000_000.0, 250.0 };
        // Supercritical
        yield return new object[] { 15_000_000.0, 500.0 };
        yield return new object[] { 20_000_000.0, 420.0 };
        // Low pressure gas
        yield return new object[] { 100_000.0, 250.0 };
    }

    [DataTestMethod]
    [DynamicData(nameof(SinglePhasePoints), DynamicDataSourceType.Method)]
    public void UpdatePH_SinglePhase_Matches_SharpFluids(double PPa, double TK)
    {
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;

        var P = Pressure.FromPascal(PPa);
        var T = Temperature.FromKelvin(TK);
        var refFluid = new Fluid(FluidList.Ammonia);
        refFluid.UpdatePT(P, T);

        double hRef = refFluid.Enthalpy!.JoulePerKilogram;
        double rhoRef = refFluid.Density!.KilogramPerCubicMeter;
        double sRef = refFluid.Entropy!.JoulePerKilogramKelvin;
        double uRef = refFluid.InternalEnergy!.JoulePerKilogram;

        var a = new AmmoniaDouble();
        a.UpdatePH(PPa, hRef);

        Assert.AreEqual(-1.0, a.Quality, "Single-phase result should report Quality=-1");
        Assert.IsTrue(Math.Abs(a.Temperature - TK) <= absTolT_K, $"T mismatch: expected={TK}, actual={a.Temperature}");
        Assert.IsTrue(RelErr(rhoRef, a.Density) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density}");
        Assert.IsTrue(RelErr(sRef, a.Entropy) <= relTol, $"s mismatch: expected={sRef}, actual={a.Entropy}");
        Assert.IsTrue(RelErr(uRef, a.InternalEnergy) <= relTol, $"u mismatch: expected={uRef}, actual={a.InternalEnergy}");
    }

    // ------------------------------------------------------------------
    // Fixed-grid two-phase: build hMix via the SharpFluids saturated endpoints (same oracle
    // pattern as UpdatePX's tests), feed it into UpdatePH, check quality/T/rho are recovered.
    // Stays within the same safe pressure range already characterized for UpdatePX (avoids
    // the 0.8-0.95 Pc ancillary-fit dip and the near-triple h/s/u sensitivity).
    // ------------------------------------------------------------------
    public static IEnumerable<object[]> TwoPhasePoints()
    {
        double[] pressuresBar = { 1, 5, 10, 20, 40, 60, 80 };
        double[] qualities = { 0.0, 0.2, 0.5, 0.8, 1.0 };

        foreach (var p in pressuresBar)
            foreach (var q in qualities)
                yield return new object[] { p, q };
    }

    [DataTestMethod]
    [DynamicData(nameof(TwoPhasePoints), DynamicDataSourceType.Method)]
    public void UpdatePH_TwoPhase_Matches_SharpFluids(double pBar, double q)
    {
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;
        const double absTolQ = 1e-4;

        var P = Pressure.FromBar(pBar);
        var L = new Fluid(FluidList.Ammonia);
        L.UpdatePX(P, 0.0);
        var V = new Fluid(FluidList.Ammonia);
        V.UpdatePX(P, 1.0);

        double rhoL = L.Density!.KilogramPerCubicMeter;
        double rhoV = V.Density!.KilogramPerCubicMeter;
        double vMixRef = (1.0 - q) / rhoL + q / rhoV;
        double rhoRef = 1.0 / vMixRef;
        double hMix = (1.0 - q) * L.Enthalpy!.JoulePerKilogram + q * V.Enthalpy!.JoulePerKilogram;

        var satFluid = new Fluid(FluidList.Ammonia);
        satFluid.Pressure = P;
        double TsatRef = satFluid.Tsat!.Kelvin;

        var a = new AmmoniaDouble();
        a.UpdatePH(P.Pascal, hMix);

        Assert.AreEqual(EngineeringFluids.Helmholtz.Phase.Phases.Twophase, a.Phase, "Phase should be Twophase");
        Assert.IsTrue(Math.Abs(a.Quality - q) <= absTolQ, $"Quality mismatch: expected={q}, actual={a.Quality}");
        Assert.IsTrue(Math.Abs(a.Temperature - TsatRef) <= absTolT_K, $"T mismatch: expected={TsatRef}, actual={a.Temperature}");
        Assert.IsTrue(RelErr(rhoRef, a.Density) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density}");
    }

    [TestMethod]
    public void UpdatePH_Throws_ForNonPositivePressure()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePH(0.0, 1_000_000.0));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePH(-1.0, 1_000_000.0));
    }

    [TestMethod]
    public void UpdatePH_Throws_ForNonFiniteEnthalpy()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePH(1_000_000.0, double.NaN));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePH(1_000_000.0, double.PositiveInfinity));
    }

    // ------------------------------------------------------------------
    // Broad, edge-biased sweep: random (P,T) single-phase states (gas/liquid/supercritical)
    // AND random (P,q) two-phase states, both converted to an enthalpy target the same way
    // a real caller would arrive at one, then verified end-to-end through UpdatePH. Learns
    // from the UpdatePT/UpdatePX sweeps: non-throwing oracle lookups (FailState checked only
    // AFTER touching lazily-evaluated properties), and a hard attempt cap so a high skip rate
    // can never hang the run.
    // ------------------------------------------------------------------
    [TestMethod]
    [TestCategory("LongRunning")]
    public void Sweep_UpdatePH_vs_SharpFluids()
    {
        const int targetPoints = 5_000;
        const int maxAttempts = targetPoints * 20;
        const double relTol = 1e-5;
        const double absTolT_K = 5e-2;

        var rng = new Random(135792468);
        int pass = 0, fail = 0, skip = 0, attempts = 0;
        var worst = new List<(double eRho, double eH, double eS, double eU, double eT, string region, double P, double q_or_T)>();

        // Stay within the ranges already validated as reliable for the underlying ancillaries
        // (see UpdatePX.cs remarks): away from the 0.8-0.95 Pc dip and the sub-1K-from-Tc
        // breakdown. The two-phase lower bound is much higher here than UpdatePX's own
        // (Ptriple*4) needs to be: UpdatePH inverts h back to a quality, and at low absolute
        // pressure rhoV is tiny, so d(specific volume)/dq ~ 1/rhoV explodes near q~0/q~1 - any
        // residual mismatch between our EOS's hL/hV and CoolProp's true values (however tiny)
        // gets amplified into real density error there. UpdatePX never hits this since it's
        // given quality directly rather than inverting it from h. First-run sweep data showed
        // this clearing up around 400 kPa; 1 MPa gives comfortable margin.
        double pLowTwoPhase = 1_000_000.0;
        double pHighTwoPhase = Pc * 0.75;

        int idx = 0;
        while (idx < targetPoints && attempts < maxAttempts)
        {
            attempts++;
            bool twoPhase = rng.NextDouble() < 0.4;

            if (twoPhase)
            {
                double PPa = Math.Exp(Math.Log(pLowTwoPhase) + rng.NextDouble() * (Math.Log(pHighTwoPhase) - Math.Log(pLowTwoPhase)));
                double q = rng.NextDouble();
                var P = Pressure.FromPascal(PPa);

                var L = new Fluid(FluidList.Ammonia);
                L.UpdatePX(P, 0.0);
                var V = new Fluid(FluidList.Ammonia);
                V.UpdatePX(P, 1.0);
                if (L.FailState || V.FailState || L.Density == null || V.Density == null || L.Enthalpy == null || V.Enthalpy == null)
                {
                    skip++;
                    continue;
                }

                var satFluid = new Fluid(FluidList.Ammonia);
                satFluid.Pressure = P;
                var Tsat = satFluid.Tsat; // lazily evaluated - read before checking FailState
                if (satFluid.FailState || Tsat == null)
                {
                    skip++;
                    continue;
                }

                double rhoL = L.Density.KilogramPerCubicMeter;
                double rhoV = V.Density.KilogramPerCubicMeter;
                double rhoRef = 1.0 / ((1.0 - q) / rhoL + q / rhoV);
                double hRef = (1.0 - q) * L.Enthalpy.JoulePerKilogram + q * V.Enthalpy.JoulePerKilogram;

                var a = new AmmoniaDouble();
                try { a.UpdatePH(PPa, hRef); }
                catch (Exception)
                {
                    fail++; idx++;
                    if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, "2P", PPa, q));
                    continue;
                }

                double eRho = RelErr(rhoRef, a.Density);
                double eQ = Math.Abs(a.Quality - q);
                double eT = Math.Abs(a.Temperature - Tsat.Kelvin);
                bool ok = eRho <= relTol && eQ <= 1e-4 && eT <= absTolT_K && a.Phase == EngineeringFluids.Helmholtz.Phase.Phases.Twophase;

                if (ok) pass++;
                else
                {
                    fail++;
                    if (worst.Count < 50) worst.Add((eRho, eQ, 0, 0, eT, "2P", PPa, q));
                }
                idx++;
            }
            else
            {
                // Single-phase: sample across gas/liquid/supercritical using the same overall
                // domain the UpdatePT correctness sweep already validated.
                double PPa = Math.Exp(Math.Log(20_000.0) + rng.NextDouble() * (Math.Log(3.0 * Pc) - Math.Log(20_000.0)));
                double TK = 195.6 + rng.NextDouble() * (700.0 - 195.6);

                var P = Pressure.FromPascal(PPa);
                var T = Temperature.FromKelvin(TK);
                var refFluid = new Fluid(FluidList.Ammonia);
                refFluid.UpdatePT(P, T);
                if (refFluid.FailState || refFluid.Density == null || refFluid.Enthalpy == null
                    || refFluid.Entropy == null || refFluid.InternalEnergy == null)
                {
                    skip++;
                    continue;
                }
                // Skip oracle points that are themselves two-phase (ambiguous single-phase target).
                if (refFluid.Quality >= 0.0 && refFluid.Quality <= 1.0)
                {
                    skip++;
                    continue;
                }

                double hRef = refFluid.Enthalpy.JoulePerKilogram;
                double rhoRef = refFluid.Density.KilogramPerCubicMeter;
                double sRef = refFluid.Entropy.JoulePerKilogramKelvin;
                double uRef = refFluid.InternalEnergy.JoulePerKilogram;

                var a = new AmmoniaDouble();
                try { a.UpdatePH(PPa, hRef); }
                catch (Exception)
                {
                    fail++; idx++;
                    if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, "1P", PPa, TK));
                    continue;
                }

                double eRho = RelErr(rhoRef, a.Density);
                double eS = RelErr(sRef, a.Entropy);
                double eU = RelErr(uRef, a.InternalEnergy);
                double eT = Math.Abs(a.Temperature - TK);
                bool ok = eRho <= relTol && eS <= relTol && eU <= relTol && eT <= absTolT_K;

                if (ok) pass++;
                else
                {
                    fail++;
                    if (worst.Count < 50) worst.Add((eRho, 0, eS, eU, eT, "1P", PPa, TK));
                }
                idx++;
            }
        }

        string summary = $"UpdatePH sweep: attempts={attempts}, tested={idx}, pass={pass}, fail={fail}, skip={skip}";
        Assert.IsTrue(idx >= targetPoints, $"Hit the attempt cap before reaching {targetPoints} tested points ({summary}) - skip rate is unexpectedly high.");

        if (fail > 0)
        {
            double NormMax((double eRho, double eH, double eS, double eU, double eT, string region, double P, double q_or_T) w)
                => Math.Max(Math.Max(w.eRho, w.eH), Math.Max(w.eS, w.eU)) / relTol;

            worst.Sort((x, y) => NormMax(y).CompareTo(NormMax(x)));
            summary += "\nWorst points:\n";
            foreach (var w in worst)
                summary += $"  [{w.region}] P={w.P:G6} Pa, q/T={w.q_or_T:G6}, eRho={w.eRho:G3}, e2={w.eH:G3}, eS={w.eS:G3}, eU={w.eU:G3}, eT={w.eT:G3}\n";
            Assert.Fail(summary);
        }
    }
}
