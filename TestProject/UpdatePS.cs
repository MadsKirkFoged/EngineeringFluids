using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;

namespace TestProject;

[TestClass]
public class CoolPropOracle_UpdatePS_Tests
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
    // entropy back into our UpdatePS, and check it recovers the same T/rho/h/u.
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
    public void UpdatePS_SinglePhase_Matches_SharpFluids(double PPa, double TK)
    {
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;

        var P = EngineeringUnits.Pressure.FromPascal(PPa);
        var T = EngineeringUnits.Temperature.FromKelvin(TK);
        var refFluid = new Fluid(FluidList.Ammonia);
        refFluid.UpdatePT(P, T);

        double sRef = refFluid.Entropy!.JoulePerKilogramKelvin;
        double rhoRef = refFluid.Density!.KilogramPerCubicMeter;
        double hRef = refFluid.Enthalpy!.JoulePerKilogram;
        double uRef = refFluid.InternalEnergy!.JoulePerKilogram;

        var a = new Ammonia();
        a.UpdatePS(Pressure.FromPascal(PPa), SpecificEntropy.FromJoulePerKilogramKelvin(sRef));

        Assert.AreEqual(-1.0, a.Quality, "Single-phase result should report Quality=-1");
        Assert.IsTrue(Math.Abs(a.Temperature.Kelvin - TK) <= absTolT_K, $"T mismatch: expected={TK}, actual={a.Temperature.Kelvin}");
        Assert.IsTrue(RelErr(rhoRef, a.Density.KilogramPerCubicMeter) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density.KilogramPerCubicMeter}");
        Assert.IsTrue(RelErr(hRef, a.Enthalpy.JoulePerKilogram) <= relTol, $"h mismatch: expected={hRef}, actual={a.Enthalpy.JoulePerKilogram}");
        Assert.IsTrue(RelErr(uRef, a.InternalEnergy.JoulePerKilogram) <= relTol, $"u mismatch: expected={uRef}, actual={a.InternalEnergy.JoulePerKilogram}");
    }

    // ------------------------------------------------------------------
    // Fixed-grid two-phase: build sMix via the SharpFluids saturated endpoints (same oracle
    // pattern as UpdatePH's tests), feed it into UpdatePS, check quality/T/rho are recovered.
    // Stays within the same safe pressure range already characterized for UpdatePX/UpdatePH.
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
    public void UpdatePS_TwoPhase_Matches_SharpFluids(double pBar, double q)
    {
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;
        const double absTolQ = 1e-4;

        var P = EngineeringUnits.Pressure.FromBar(pBar);
        var L = new Fluid(FluidList.Ammonia);
        L.UpdatePX(P, 0.0);
        var V = new Fluid(FluidList.Ammonia);
        V.UpdatePX(P, 1.0);

        double rhoL = L.Density!.KilogramPerCubicMeter;
        double rhoV = V.Density!.KilogramPerCubicMeter;
        double vMixRef = (1.0 - q) / rhoL + q / rhoV;
        double rhoRef = 1.0 / vMixRef;
        double sMix = (1.0 - q) * L.Entropy!.JoulePerKilogramKelvin + q * V.Entropy!.JoulePerKilogramKelvin;

        var satFluid = new Fluid(FluidList.Ammonia);
        satFluid.Pressure = P;
        double TsatRef = satFluid.Tsat!.Kelvin;

        var a = new Ammonia();
        a.UpdatePS(Pressure.FromPascal(P.Pascal), SpecificEntropy.FromJoulePerKilogramKelvin(sMix));

        Assert.AreEqual(EngineeringFluids.Helmholtz.Phase.Phases.Twophase, a.Phase, "Phase should be Twophase");
        Assert.IsTrue(Math.Abs(a.Quality - q) <= absTolQ, $"Quality mismatch: expected={q}, actual={a.Quality}");
        Assert.IsTrue(Math.Abs(a.Temperature.Kelvin - TsatRef) <= absTolT_K, $"T mismatch: expected={TsatRef}, actual={a.Temperature.Kelvin}");
        Assert.IsTrue(RelErr(rhoRef, a.Density.KilogramPerCubicMeter) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density.KilogramPerCubicMeter}");
    }

    [TestMethod]
    public void UpdatePS_Throws_ForNonPositivePressure()
    {
        var a = new Ammonia();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePS(Pressure.FromPascal(0.0), SpecificEntropy.FromJoulePerKilogramKelvin(1000.0)));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePS(Pressure.FromPascal(-1.0), SpecificEntropy.FromJoulePerKilogramKelvin(1000.0)));
    }

    [TestMethod]
    public void UpdatePS_Throws_ForNonFiniteEntropy()
    {
        var a = new Ammonia();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePS(Pressure.FromPascal(1_000_000.0), SpecificEntropy.FromJoulePerKilogramKelvin(double.NaN)));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePS(Pressure.FromPascal(1_000_000.0), SpecificEntropy.FromJoulePerKilogramKelvin(double.PositiveInfinity)));
    }

    // ------------------------------------------------------------------
    // Broad, edge-biased sweep: random (P,T) single-phase states AND random (P,q) two-phase
    // states, both converted to an entropy target the same way a real caller would arrive at
    // one, then verified end-to-end through UpdatePS. Same non-throwing oracle pattern and
    // attempt cap as the UpdatePH sweep.
    // ------------------------------------------------------------------
    [TestMethod]
    [TestCategory("LongRunning")]
    public void Sweep_UpdatePS_vs_SharpFluids()
    {
        const int targetPoints = 5_000;
        const int maxAttempts = targetPoints * 20;
        const double relTol = 1e-5;
        const double absTolT_K = 5e-2;

        var rng = new Random(246813579);
        int pass = 0, fail = 0, skip = 0, attempts = 0;
        var worst = new List<(double eRho, double eH, double eS, double eU, double eT, string region, double P, double q_or_T)>();

        // Same reliable ranges already validated for UpdatePX/UpdatePH (see their remarks):
        // away from the 0.8-0.95 Pc ancillary-fit dip and the sub-1K-from-Tc breakdown, and a
        // two-phase lower bound high enough to stay clear of the low-pressure amplification
        // that comes from inverting a mixture property back to a quality (1/rhoV blows up as
        // rhoV -> 0 near q~0/q~1 at low absolute pressure) - inherent to UpdatePS the same way
        // it is to UpdatePH, since both invert a mixture property rather than being given
        // quality directly the way UpdatePX is.
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
                var P = EngineeringUnits.Pressure.FromPascal(PPa);

                var L = new Fluid(FluidList.Ammonia);
                L.UpdatePX(P, 0.0);
                var V = new Fluid(FluidList.Ammonia);
                V.UpdatePX(P, 1.0);
                if (L.FailState || V.FailState || L.Density == null || V.Density == null || L.Entropy == null || V.Entropy == null)
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
                double sRef = (1.0 - q) * L.Entropy.JoulePerKilogramKelvin + q * V.Entropy.JoulePerKilogramKelvin;

                var a = new Ammonia();
                try { a.UpdatePS(Pressure.FromPascal(PPa), SpecificEntropy.FromJoulePerKilogramKelvin(sRef)); }
                catch (Exception)
                {
                    fail++; idx++;
                    if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, "2P", PPa, q));
                    continue;
                }

                double eRho = RelErr(rhoRef, a.Density.KilogramPerCubicMeter);
                double eQ = Math.Abs(a.Quality - q);
                double eT = Math.Abs(a.Temperature.Kelvin - Tsat.Kelvin);
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
                double PPa = Math.Exp(Math.Log(20_000.0) + rng.NextDouble() * (Math.Log(3.0 * Pc) - Math.Log(20_000.0)));
                double TK = 195.6 + rng.NextDouble() * (700.0 - 195.6);

                // Within ~1K of Tc, entropy is dominated by the ideal-gas term (delta is tiny
                // there for the low/moderate pressures this sweep covers) but the two EOS
                // implementations' near-critical residual contributions diverge from each
                // other by a small, genuine amount - confirmed by direct reproduction: fed its
                // OWN entropy back in, this EOS recovers T to 1.5e-7 K, so the root-find itself
                // is exact; the mismatch is purely this-EOS-vs-CoolProp disagreement right at
                // Tc, the same known reduced-reliability band already documented for
                // UpdatePX/UpdatePH (see their remarks) applied to entropy here too.
                if (Math.Abs(TK - Tc) < 1.0)
                {
                    skip++;
                    continue;
                }

                var P = EngineeringUnits.Pressure.FromPascal(PPa);
                var T = EngineeringUnits.Temperature.FromKelvin(TK);
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

                double sRef = refFluid.Entropy.JoulePerKilogramKelvin;
                double rhoRef = refFluid.Density.KilogramPerCubicMeter;
                double hRef = refFluid.Enthalpy.JoulePerKilogram;
                double uRef = refFluid.InternalEnergy.JoulePerKilogram;

                var a = new Ammonia();
                try { a.UpdatePS(Pressure.FromPascal(PPa), SpecificEntropy.FromJoulePerKilogramKelvin(sRef)); }
                catch (Exception)
                {
                    fail++; idx++;
                    if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, "1P", PPa, TK));
                    continue;
                }

                double eRho = RelErr(rhoRef, a.Density.KilogramPerCubicMeter);
                double eH = RelErr(hRef, a.Enthalpy.JoulePerKilogram);
                double eU = RelErr(uRef, a.InternalEnergy.JoulePerKilogram);
                double eT = Math.Abs(a.Temperature.Kelvin - TK);
                bool ok = eRho <= relTol && eH <= relTol && eU <= relTol && eT <= absTolT_K;

                if (ok) pass++;
                else
                {
                    fail++;
                    if (worst.Count < 50) worst.Add((eRho, eH, 0, eU, eT, "1P", PPa, TK));
                }
                idx++;
            }
        }

        string summary = $"UpdatePS sweep: attempts={attempts}, tested={idx}, pass={pass}, fail={fail}, skip={skip}";
        Assert.IsTrue(idx >= targetPoints, $"Hit the attempt cap before reaching {targetPoints} tested points ({summary}) - skip rate is unexpectedly high.");

        if (fail > 0)
        {
            double NormMax((double eRho, double eH, double eS, double eU, double eT, string region, double P, double q_or_T) w)
                => Math.Max(Math.Max(w.eRho, w.eH), Math.Max(w.eS, w.eU)) / relTol;

            worst.Sort((x, y) => NormMax(y).CompareTo(NormMax(x)));
            summary += "\nWorst points:\n";
            foreach (var w in worst)
                summary += $"  [{w.region}] P={w.P:G6} Pa, q/T={w.q_or_T:G6}, eRho={w.eRho:G3}, eH={w.eH:G3}, eU={w.eU:G3}, eT={w.eT:G3}\n";
            Assert.Fail(summary);
        }
    }
}
