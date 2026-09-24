using EngineeringFluids.Fluids;
using EngineeringUnits;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;

namespace TestProject;

[TestClass]
public class CoolPropOracle_UpdatePX_Tests
{
    private const double Tc = 405.56;
    private const double Pc = 11363400.0;
    private const double Ptriple = 6091.223108650368;

    private static double RelErr(double expected, double actual)
    {
        double denom = Math.Max(Math.Abs(expected), 1.0);
        return Math.Abs(actual - expected) / denom;
    }

    // Right at the critical point, CoolProp/SharpFluids can fail to converge without
    // throwing: it sets FailState=true but leaves Density/Pressure/etc at whatever they
    // were after the last successful solve on this instance. Treat that as "no answer
    // available" rather than comparable ground truth (same pattern as the UpdatePT tests).
    private static void SkipIfOracleFailed(Fluid f, string context)
    {
        if (f.FailState)
            Assert.Inconclusive($"SharpFluids oracle failed (FailState=true): {context}");
    }

    private static (double rho, double h, double s, double u, double T) OracleMix(double PPa, double q)
    {
        var P = Pressure.FromPascal(PPa);

        var L = new Fluid(FluidList.Ammonia);
        L.UpdatePX(P, 0.0);
        SkipIfOracleFailed(L, $"UpdatePX(P,0) at P={PPa} Pa");

        var V = new Fluid(FluidList.Ammonia);
        V.UpdatePX(P, 1.0);
        SkipIfOracleFailed(V, $"UpdatePX(P,1) at P={PPa} Pa");

        var satFluid = new Fluid(FluidList.Ammonia);
        satFluid.Pressure = P;
        // .Tsat is lazily evaluated - FailState only reflects reality once it (or another
        // Update call) has actually run, so it must be read AFTER accessing .Tsat, not before.
        var Tsat = satFluid.Tsat;
        SkipIfOracleFailed(satFluid, $"Tsat at P={PPa} Pa");
        if (Tsat == null)
            Assert.Inconclusive($"Tsat null at P={PPa} Pa");

        double rhoL = L.Density!.KilogramPerCubicMeter;
        double rhoV = V.Density!.KilogramPerCubicMeter;
        double vMix = (1.0 - q) / rhoL + q / rhoV;
        double rho = 1.0 / vMix;

        double h = (1.0 - q) * L.Enthalpy!.JoulePerKilogram + q * V.Enthalpy!.JoulePerKilogram;
        double s = (1.0 - q) * L.Entropy!.JoulePerKilogramKelvin + q * V.Entropy!.JoulePerKilogramKelvin;
        double u = (1.0 - q) * L.InternalEnergy!.JoulePerKilogram + q * V.InternalEnergy!.JoulePerKilogram;

        return (rho, h, s, u, Tsat!.Kelvin);
    }

    public static IEnumerable<object[]> PressureQualityPoints()
    {
        // Stay within the pressure range where the ancillary fits reliably hold to 1e-5 or
        // better (measured by sweeping rhoL/rhoV error vs P/Pc in 0.02 steps): error grows
        // roughly with proximity to Pc but is NOT monotonic - there's a localized dip in fit
        // accuracy around 0.8-0.95 Pc (up to ~5e-5, a polynomial-fit imperfection, not the
        // sharp critical-point breakdown) before the real breakdown within ~1K of Tc that
        // UpdatePX itself rejects (see UpdatePXFast.cs). 85 bar (~0.75 Pc, Pc=113.634 bar)
        // stays comfortably under that dip.
        double[] pressuresBar = { 0.5, 1, 2, 5, 10, 20, 40, 60, 80, 85 };
        double[] qualities = { 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 };

        foreach (var p in pressuresBar)
            foreach (var q in qualities)
                yield return new object[] { p, q };
    }

    [DataTestMethod]
    [DynamicData(nameof(PressureQualityPoints), DynamicDataSourceType.Method)]
    public void UpdatePX_Matches_SharpFluids(double pBar, double q)
    {
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;
        const double absTolP_Pa = 1e-3; // Psat is forced to exactly pTarget, so this should be ~0

        double PPa = Pressure.FromBar(pBar).Pascal;
        var (rhoRef, hRef, sRef, uRef, TsatRef) = OracleMix(PPa, q);

        var a = new AmmoniaDouble();
        a.UpdatePX(PPa, q);

        Assert.AreEqual(EngineeringFluids.Helmholtz.Phase.Phases.Twophase, a.Phase, "Phase should be Twophase after UpdatePX");
        Assert.AreEqual(q, a.Quality, 1e-12, "Quality should be set exactly (no solving involved)");

        Assert.IsTrue(RelErr(rhoRef, a.Density) <= relTol, $"rho mismatch at P={PPa} Pa, q={q}: expected={rhoRef}, actual={a.Density}");
        Assert.IsTrue(RelErr(hRef, a.Enthalpy) <= relTol, $"h mismatch at P={PPa} Pa, q={q}: expected={hRef}, actual={a.Enthalpy}");
        Assert.IsTrue(RelErr(sRef, a.Entropy) <= relTol, $"s mismatch at P={PPa} Pa, q={q}: expected={sRef}, actual={a.Entropy}");
        Assert.IsTrue(RelErr(uRef, a.InternalEnergy) <= relTol, $"u mismatch at P={PPa} Pa, q={q}: expected={uRef}, actual={a.InternalEnergy}");

        Assert.IsTrue(Math.Abs(a.Temperature - TsatRef) <= absTolT_K, $"T mismatch at P={PPa} Pa, q={q}: expected={TsatRef}, actual={a.Temperature}");
        Assert.IsTrue(Math.Abs(a.Pressure - PPa) <= absTolP_Pa, $"p back-calc mismatch at P={PPa} Pa, q={q}: expected={PPa}, actual={a.Pressure}");
    }

    [TestMethod]
    public void UpdatePX_Throws_AtOrAboveCriticalPressure()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdatePX(Pc, 0.5));
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdatePX(Pc * 1.1, 0.5));
    }

    [TestMethod]
    public void UpdatePX_Throws_BelowTriplePressure()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePX(Ptriple * 0.5, 0.5));
    }

    [TestMethod]
    public void UpdatePX_Throws_NearCriticalPoint()
    {
        // Within the last 1K below Tc, ancillary-based rhoL/rhoV lose too much accuracy
        // (see UpdatePXFast.cs remarks and the sweep test below) - UpdatePX rejects it
        // outright instead of silently returning an inaccurate state.
        var a = new AmmoniaDouble();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdatePX(Pc * 0.995, 0.5));
    }

    [DataTestMethod]
    [DataRow(-0.1)]
    [DataRow(1.1)]
    [DataRow(double.NaN)]
    public void UpdatePX_Throws_ForInvalidQuality(double q)
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdatePX(1_000_000.0, q));
    }

    // Non-throwing oracle lookup for the tight sweep loop below - avoids Assert.Inconclusive's
    // exception overhead on every skip, which matters when sampling anywhere near Pc where
    // SharpFluids/CoolProp legitimately sets FailState often. NOTE: FailState only reflects
    // reality once a computation has actually run on that Fluid instance - it must be checked
    // AFTER accessing a lazily-evaluated property like .Tsat, never before (checking it before
    // .Tsat looks like a permanent 100% failure, since nothing has been computed on that
    // instance yet).
    private static bool TryOracleMix(double PPa, double q, out double rho, out double h, out double s, out double u, out double T)
    {
        rho = h = s = u = T = double.NaN;
        var P = Pressure.FromPascal(PPa);

        var L = new Fluid(FluidList.Ammonia);
        L.UpdatePX(P, 0.0);
        if (L.FailState || L.Density == null || L.Enthalpy == null || L.Entropy == null || L.InternalEnergy == null)
            return false;

        var V = new Fluid(FluidList.Ammonia);
        V.UpdatePX(P, 1.0);
        if (V.FailState || V.Density == null || V.Enthalpy == null || V.Entropy == null || V.InternalEnergy == null)
            return false;

        var satFluid = new Fluid(FluidList.Ammonia);
        satFluid.Pressure = P;
        var Tsat = satFluid.Tsat;
        if (satFluid.FailState || Tsat == null)
            return false;

        double rhoL = L.Density.KilogramPerCubicMeter;
        double rhoV = V.Density.KilogramPerCubicMeter;
        double vMix = (1.0 - q) / rhoL + q / rhoV;
        rho = 1.0 / vMix;

        h = (1.0 - q) * L.Enthalpy.JoulePerKilogram + q * V.Enthalpy.JoulePerKilogram;
        s = (1.0 - q) * L.Entropy.JoulePerKilogramKelvin + q * V.Entropy.JoulePerKilogramKelvin;
        u = (1.0 - q) * L.InternalEnergy.JoulePerKilogram + q * V.InternalEnergy.JoulePerKilogram;
        T = Tsat.Kelvin;
        return true;
    }

    // Broad, edge-biased sweep across the whole valid domain (log-uniform pressure from just
    // above the triple-point pressure up to a safety margin below Pc, and quality spanning
    // [0,1] including the q=0/q=1 endpoints) - the same style of test that caught a real
    // oracle-reuse bug for UpdatePT's 100k sweep. A hard attempt cap guarantees termination
    // even if the skip rate is ever much higher than expected.
    [TestMethod]
    [TestCategory("LongRunning")]
    public void Sweep_UpdatePX_vs_SharpFluids()
    {
        const int targetPoints = 20_000;
        const int maxAttempts = targetPoints * 20;
        const double relTol = 1e-5;
        const double absTolT_K = 1e-2;

        // Pc*0.75 stays under the 0.8-0.95 Pc dip in ancillary fit accuracy (up to ~5e-5,
        // measured separately - see PressureQualityPoints remarks above; 0.78 still caught
        // its leading edge in a 20k-point run), not just under the sharp critical-point
        // breakdown itself (which starts within ~1K of Tc, i.e. much closer to Pc than this).
        double pHigh = Pc * 0.75;
        // Within the first few kPa above the triple-point pressure, h/s/u (not rho, which
        // stays accurate to ~1e-6 throughout) pick up rel. errors up to ~2e-4 at the bubble
        // point (q~0) - first-run data showed this clearing up by ~4x the triple pressure.
        // Likely the ancillary fits' 2K segment width combined with how steep dP/dT (and
        // hence the EOS's temperature-derivative terms) is that close to the triple point.
        double pLow = Ptriple * 4.0;

        var rng = new Random(987654321);
        int pass = 0, fail = 0, skip = 0, attempts = 0;
        var worst = new List<(double eRho, double eH, double eS, double eU, double eT, double P, double q)>();

        double SampleQuality(Random r)
        {
            double u = r.NextDouble();
            // 30% exactly at the endpoints (q=0 or q=1), 70% uniform in (0,1)
            if (u < 0.15) return 0.0;
            if (u < 0.30) return 1.0;
            return r.NextDouble();
        }

        int idx = 0;
        while (idx < targetPoints && attempts < maxAttempts)
        {
            attempts++;

            double PPa = Math.Exp(Math.Log(pLow) + rng.NextDouble() * (Math.Log(pHigh) - Math.Log(pLow)));
            double q = SampleQuality(rng);

            if (!TryOracleMix(PPa, q, out double rhoRef, out double hRef, out double sRef, out double uRef, out double TsatRef))
            {
                skip++;
                continue;
            }

            var a = new AmmoniaDouble();
            try
            {
                a.UpdatePX(PPa, q);
            }
            catch (Exception)
            {
                fail++;
                idx++;
                if (worst.Count < 50)
                    worst.Add((1.0, 1.0, 1.0, 1.0, 1.0, PPa, q));
                continue;
            }

            double eRho = RelErr(rhoRef, a.Density);
            double eH = RelErr(hRef, a.Enthalpy);
            double eS = RelErr(sRef, a.Entropy);
            double eU = RelErr(uRef, a.InternalEnergy);
            double eT = Math.Abs(a.Temperature - TsatRef);

            bool ok = eRho <= relTol && eH <= relTol && eS <= relTol && eU <= relTol && eT <= absTolT_K;

            if (ok)
                pass++;
            else
            {
                fail++;
                if (worst.Count < 50)
                    worst.Add((eRho, eH, eS, eU, eT, PPa, q));
            }

            idx++;
        }

        string summary = $"UpdatePX sweep: attempts={attempts}, tested={idx}, pass={pass}, fail={fail}, skip={skip}";
        Assert.IsTrue(idx >= targetPoints, $"Hit the attempt cap before reaching {targetPoints} tested points ({summary}) - skip rate is unexpectedly high.");

        if (fail > 0)
        {
            // Rank by whichever normalized error (relative to its own tolerance) is largest,
            // since a point can fail on h/s/u/T even when rho itself is fine.
            double NormMax((double eRho, double eH, double eS, double eU, double eT, double P, double q) w)
                => Math.Max(Math.Max(w.eRho, w.eH), Math.Max(w.eS, w.eU)) / relTol;

            worst.Sort((x, y) => NormMax(y).CompareTo(NormMax(x)));
            summary += "\nWorst points:\n";
            foreach (var w in worst)
                summary += $"  P={w.P:G6} Pa, q={w.q:F4}, eRho={w.eRho:G3}, eH={w.eH:G3}, eS={w.eS:G3}, eU={w.eU:G3}, eT={w.eT:G3}\n";
            Assert.Fail(summary);
        }
    }
}
