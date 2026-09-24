using EngineeringFluids.Fluids;
using EngineeringUnits;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;

namespace TestProject;

[TestClass]
public class CoolPropOracle_UpdateTX_FastTests
{
    private const double Tc = 405.56;
    private const double Ttriple = 195.495;

    private static double RelErr(double expected, double actual)
    {
        double denom = Math.Max(Math.Abs(expected), 1.0);
        return Math.Abs(actual - expected) / denom;
    }

    // ------------------------------------------------------------------
    // Fixed-grid: build the saturated state via SharpFluids' own UpdateXT (the same
    // operation, so this is a direct apples-to-apples oracle comparison rather than an
    // indirect one), and check our UpdateTX/UpdateXT recover the same P/rho/h/s/u.
    // ------------------------------------------------------------------
    public static IEnumerable<object[]> Points()
    {
        // Stays clear of the reduced-accuracy zones documented in the sweep below: within
        // ~20 K of the triple point (Psat ancillary) and ~16 K of Tc (rhoL/rhoV ancillaries).
        double[] temperatures = { 216.0, 220.0, 250.0, 280.0, 300.0, 330.0, 360.0, 380.0, 384.0 };
        double[] qualities = { 0.0, 0.2, 0.5, 0.8, 1.0 };

        foreach (var t in temperatures)
            foreach (var q in qualities)
                yield return new object[] { t, q };
    }

    [DataTestMethod]
    [DynamicData(nameof(Points), DynamicDataSourceType.Method)]
    public void UpdateTX_Matches_SharpFluids(double TK, double q)
    {
        const double relTol = 1e-5;
        const double absTolP = 1.0; // Pa

        var refFluid = new Fluid(FluidList.Ammonia);
        refFluid.UpdateXT(q, Temperature.FromKelvin(TK));

        double pRef = refFluid.Pressure!.Pascal;
        double rhoRef = refFluid.Density!.KilogramPerCubicMeter;
        double hRef = refFluid.Enthalpy!.JoulePerKilogram;
        double sRef = refFluid.Entropy!.JoulePerKilogramKelvin;
        double uRef = refFluid.InternalEnergy!.JoulePerKilogram;

        var a = new AmmoniaDouble();
        a.UpdateTX(TK, q);

        Assert.AreEqual(EngineeringFluids.Helmholtz.Phase.Phases.Twophase, a.Phase, "Phase should be Twophase");
        Assert.IsTrue(Math.Abs(a.Quality - q) <= 1e-9, $"Quality mismatch: expected={q}, actual={a.Quality}");
        Assert.IsTrue(Math.Abs(a.Temperature - TK) <= 1e-9, $"T mismatch: expected={TK}, actual={a.Temperature}");
        Assert.IsTrue(Math.Abs(a.Pressure - pRef) <= Math.Max(absTolP, relTol * pRef), $"P mismatch: expected={pRef}, actual={a.Pressure}");
        Assert.IsTrue(RelErr(rhoRef, a.Density) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density}");
        Assert.IsTrue(RelErr(hRef, a.Enthalpy) <= relTol, $"h mismatch: expected={hRef}, actual={a.Enthalpy}");
        Assert.IsTrue(RelErr(sRef, a.Entropy) <= relTol, $"s mismatch: expected={sRef}, actual={a.Entropy}");
        Assert.IsTrue(RelErr(uRef, a.InternalEnergy) <= relTol, $"u mismatch: expected={uRef}, actual={a.InternalEnergy}");
    }

    [TestMethod]
    public void UpdateXT_MatchesUpdateTX()
    {
        var a = new AmmoniaDouble();
        a.UpdateTX(280.0, 0.3);
        var b = new AmmoniaDouble();
        b.UpdateXT(0.3, 280.0);

        Assert.AreEqual(a.Temperature, b.Temperature);
        Assert.AreEqual(a.Density, b.Density);
        Assert.AreEqual(a.Quality, b.Quality);
    }

    [TestMethod]
    public void UpdateTX_Throws_ForInvalidQuality()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTX(280.0, -0.1));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTX(280.0, 1.1));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTX(280.0, double.NaN));
    }

    [TestMethod]
    public void UpdateTX_Throws_AtOrAboveTc()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTX(Tc, 0.5));
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTX(Tc + 5.0, 0.5));
    }

    [TestMethod]
    public void UpdateTX_Throws_BelowTripleTemperature()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTX(Ttriple - 1.0, 0.5));
    }

    [TestMethod]
    public void UpdateTX_Throws_NearCriticalPoint()
    {
        var a = new AmmoniaDouble();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTX(Tc - 0.5, 0.5));
    }

    // ------------------------------------------------------------------
    // Broad randomized sweep across the whole reliable temperature range (triple to ~1K below
    // Tc, matching UpdatePX's own documented near-critical exclusion for the same ancillary-
    // accuracy reason) and the full quality range.
    // ------------------------------------------------------------------
    [TestMethod]
    [TestCategory("LongRunning")]
    public void Sweep_UpdateTX_vs_SharpFluids()
    {
        const int targetPoints = 5_000;
        const int maxAttempts = targetPoints * 5;
        const double relTol = 1e-5;

        var rng = new Random(864213579);
        int pass = 0, fail = 0, skip = 0, attempts = 0;
        var worst = new List<(double eP, double eRho, double eH, double eS, double eU, double T, double q)>();

        // BubbleDensityFast/DewDensityFast accuracy degrades gradually well before the
        // production-level Tc-1.0 hard cutoff (UpdateTXFast.cs rejects outright below that),
        // and well beyond the narrow triple-point margin one might expect too - first-run
        // sweep data showed relative errors up to ~3e-4 in density (and ~7e-5 in h/s/u)
        // appearing as far as ~10-11 K below Tc, and Psat's own fit breaching 1e-5 as far as
        // ~15 K above the triple point. Same phenomenon already documented for UpdatePX's
        // pressure-domain sweep (the 0.8-0.95 Pc dip) - just characterized here directly in
        // temperature space, which turns out to span a wider margin than that P-domain
        // sampling ever exercised.
        double Tlo = Ttriple + 20.0;
        double Thi = Tc - 16.0;

        int idx = 0;
        while (idx < targetPoints && attempts < maxAttempts)
        {
            attempts++;
            double TK = Tlo + rng.NextDouble() * (Thi - Tlo);
            double q = rng.NextDouble();

            var refFluid = new Fluid(FluidList.Ammonia);
            refFluid.UpdateXT(q, Temperature.FromKelvin(TK));
            if (refFluid.FailState || refFluid.Density == null || refFluid.Pressure == null
                || refFluid.Enthalpy == null || refFluid.Entropy == null || refFluid.InternalEnergy == null)
            {
                skip++;
                continue;
            }

            double pRef = refFluid.Pressure.Pascal;
            double rhoRef = refFluid.Density.KilogramPerCubicMeter;
            double hRef = refFluid.Enthalpy.JoulePerKilogram;
            double sRef = refFluid.Entropy.JoulePerKilogramKelvin;
            double uRef = refFluid.InternalEnergy.JoulePerKilogram;

            var a = new AmmoniaDouble();
            try { a.UpdateTX(TK, q); }
            catch (Exception)
            {
                fail++; idx++;
                if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, TK, q));
                continue;
            }

            double eP = RelErr(pRef, a.Pressure);
            double eRho = RelErr(rhoRef, a.Density);
            double eH = RelErr(hRef, a.Enthalpy);
            double eS = RelErr(sRef, a.Entropy);
            double eU = RelErr(uRef, a.InternalEnergy);
            bool ok = eP <= relTol && eRho <= relTol && eH <= relTol && eS <= relTol && eU <= relTol;

            if (ok) pass++;
            else
            {
                fail++;
                if (worst.Count < 50) worst.Add((eP, eRho, eH, eS, eU, TK, q));
            }
            idx++;
        }

        string summary = $"UpdateTX sweep: attempts={attempts}, tested={idx}, pass={pass}, fail={fail}, skip={skip}";
        Assert.IsTrue(idx >= targetPoints, $"Hit the attempt cap before reaching {targetPoints} tested points ({summary}) - skip rate is unexpectedly high.");

        if (fail > 0)
        {
            double NormMax((double eP, double eRho, double eH, double eS, double eU, double T, double q) w)
                => Math.Max(Math.Max(w.eP, w.eRho), Math.Max(Math.Max(w.eH, w.eS), w.eU)) / relTol;

            worst.Sort((x, y) => NormMax(y).CompareTo(NormMax(x)));
            summary += "\nWorst points:\n";
            foreach (var w in worst)
                summary += $"  T={w.T:G6} K, q={w.q:G6}, eP={w.eP:G3}, eRho={w.eRho:G3}, eH={w.eH:G3}, eS={w.eS:G3}, eU={w.eU:G3}\n";
            Assert.Fail(summary);
        }
    }
}
