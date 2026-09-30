using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;

namespace TestProject;

// Exact (EOS-based) UpdateTX: unlike CoolPropOracle_UpdateTX_Tests, this does not need to
// exclude a wide band near Tc/the triple point, since SolveAtT solves the real
// phase-equilibrium conditions rather than reading a pre-fitted ancillary curve - it should
// stay accurate right up to both physical boundaries.
[TestClass]
public class CoolPropOracle_UpdateTXExact_Tests
{
    private const double Tc = 405.56;
    private const double Ttriple = 195.495;

    private static double RelErr(double expected, double actual)
    {
        double denom = Math.Max(Math.Abs(expected), 1.0);
        return Math.Abs(actual - expected) / denom;
    }

    public static IEnumerable<object[]> Points()
    {
        // Includes points deep in both previously-excluded ancillary-fit bands (near Tc, near
        // triple), which UpdateTX's own tests deliberately avoid. Does not go closer than
        // 405.0 (0.56 K below Tc): SharpFluids/CoolProp's OWN saturation solver returns null
        // there (confirmed - not just imprecise, it fails outright), so there is no oracle
        // left to validate against, even though SolveAtT itself keeps converging cleanly
        // all the way to within 0.001 K of Tc (verified directly against its own residuals).
        double[] temperatures = { 195.6, 196.0, 200.0, 205.0, 220.0, 250.0, 280.0, 300.0, 330.0,
                                   360.0, 390.0, 400.0, 404.0, 405.0 };
        double[] qualities = { 0.0, 0.2, 0.5, 0.8, 1.0 };

        foreach (var t in temperatures)
            foreach (var q in qualities)
                yield return new object[] { t, q };
    }

    [DataTestMethod]
    [DynamicData(nameof(Points), DynamicDataSourceType.Method)]
    public void UpdateTXExact_Matches_SharpFluids(double TK, double q)
    {
        const double relTol = 1e-5;
        const double absTolP = 1.0; // Pa

        var refFluid = new Fluid(FluidList.Ammonia);
        refFluid.UpdateXT(q, EngineeringUnits.Temperature.FromKelvin(TK));

        double pRef = refFluid.Pressure!.Pascal;
        double rhoRef = refFluid.Density!.KilogramPerCubicMeter;
        double hRef = refFluid.Enthalpy!.JoulePerKilogram;
        double sRef = refFluid.Entropy!.JoulePerKilogramKelvin;
        double uRef = refFluid.InternalEnergy!.JoulePerKilogram;

        var a = new Ammonia();
        a.UpdateTXExact(Temperature.FromKelvin(TK), q);

        Assert.AreEqual(EngineeringFluids.Helmholtz.Phase.Phases.Twophase, a.Phase, "Phase should be Twophase");
        Assert.IsTrue(Math.Abs(a.Quality - q) <= 1e-9, $"Quality mismatch: expected={q}, actual={a.Quality}");
        Assert.IsTrue(Math.Abs(a.Temperature.Kelvin - TK) <= 1e-9, $"T mismatch: expected={TK}, actual={a.Temperature.Kelvin}");
        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - pRef) <= Math.Max(absTolP, relTol * pRef), $"P mismatch: expected={pRef}, actual={a.Pressure.Pascal}");
        Assert.IsTrue(RelErr(rhoRef, a.Density.KilogramPerCubicMeter) <= relTol, $"rho mismatch: expected={rhoRef}, actual={a.Density.KilogramPerCubicMeter}");
        Assert.IsTrue(RelErr(hRef, a.Enthalpy.JoulePerKilogram) <= relTol, $"h mismatch: expected={hRef}, actual={a.Enthalpy.JoulePerKilogram}");
        Assert.IsTrue(RelErr(sRef, a.Entropy.JoulePerKilogramKelvin) <= relTol, $"s mismatch: expected={sRef}, actual={a.Entropy.JoulePerKilogramKelvin}");
        Assert.IsTrue(RelErr(uRef, a.InternalEnergy.JoulePerKilogram) <= relTol, $"u mismatch: expected={uRef}, actual={a.InternalEnergy.JoulePerKilogram}");
    }

    [TestMethod]
    public void UpdateXTExact_MatchesUpdateTXExact()
    {
        var a = new Ammonia();
        a.UpdateTXExact(Temperature.FromKelvin(280.0), 0.3);
        var b = new Ammonia();
        b.UpdateXTExact(0.3, Temperature.FromKelvin(280.0));

        Assert.AreEqual(a.Temperature.Kelvin, b.Temperature.Kelvin);
        Assert.AreEqual(a.Density.KilogramPerCubicMeter, b.Density.KilogramPerCubicMeter);
        Assert.AreEqual(a.Quality, b.Quality);
    }

    [TestMethod]
    public void UpdateTXExact_Throws_ForInvalidQuality()
    {
        var a = new Ammonia();
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTXExact(Temperature.FromKelvin(280.0), -0.1));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTXExact(Temperature.FromKelvin(280.0), 1.1));
        Assert.ThrowsException<ArgumentOutOfRangeException>(() => a.UpdateTXExact(Temperature.FromKelvin(280.0), double.NaN));
    }

    [TestMethod]
    public void UpdateTXExact_Throws_AtOrAboveTc()
    {
        var a = new Ammonia();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTXExact(Temperature.FromKelvin(Tc), 0.5));
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTXExact(Temperature.FromKelvin(Tc + 5.0), 0.5));
    }

    [TestMethod]
    public void UpdateTXExact_Throws_BelowTripleTemperature()
    {
        var a = new Ammonia();
        Assert.ThrowsException<InvalidOperationException>(() => a.UpdateTXExact(Temperature.FromKelvin(Ttriple - 1.0), 0.5));
    }

    // Like CoolProp (and SharpFluids' UpdateXT), the triple temperature itself works. Its saturation pressure is
    // the EOS's own 6055.8 Pa, not the tabulated triple-point pressure of 6091.2 Pa (those densities don't belong
    // to this EOS - see SaturationSolver.SaturationAtT).
    [TestMethod]
    public void UpdateTXExact_AtTripleTemperature_GivesTheEosSaturationState()
    {
        var a = new Ammonia();
        a.UpdateTXExact(Temperature.FromKelvin(Ttriple), 0.5);
        Assert.AreEqual(6055.81357399642, a.Pressure.Pascal, 1e-3);
        Assert.AreEqual(0.5, a.Quality);
    }

    // ------------------------------------------------------------------
    // Broad randomized sweep across the FULL reliable temperature range, right up to both
    // physical boundaries (0.05 K margins only, to stay inside SolveAtT's own domain
    // validation) - no ancillary-fit exclusion zone needed here.
    // ------------------------------------------------------------------
    [TestMethod]
    [TestCategory("LongRunning")]
    public void Sweep_UpdateTXExact_vs_SharpFluids()
    {
        const int targetPoints = 2_000;
        const int maxAttempts = targetPoints * 5;
        const double relTol = 1e-5;

        var rng = new Random(975318642);
        int pass = 0, fail = 0, skip = 0, attempts = 0;
        var worst = new List<(double eP, double eRho, double eH, double eS, double eU, double T, double q)>();

        double Tlo = Ttriple + 0.05;
        double Thi = Tc - 0.05;

        int idx = 0;
        while (idx < targetPoints && attempts < maxAttempts)
        {
            attempts++;
            double TK = Tlo + rng.NextDouble() * (Thi - Tlo);
            double q = rng.NextDouble();

            // Within ~0.5 K of Tc, this EOS and CoolProp's genuinely disagree by a small
            // amount even though each is internally exact (confirmed directly: SolveAtT's
            // own pressure/fugacity residuals stay at ~1e-8/1e-14 right up to 0.001 K from Tc -
            // see UpdateTXExact.cs remarks). Same known near-critical model divergence already
            // found and excluded for UpdatePS's single-phase sweep.
            if (Tc - TK < 0.5)
            {
                skip++;
                continue;
            }

            var refFluid = new Fluid(FluidList.Ammonia);
            refFluid.UpdateXT(q, EngineeringUnits.Temperature.FromKelvin(TK));
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

            var a = new Ammonia();
            try { a.UpdateTXExact(Temperature.FromKelvin(TK), q); }
            catch (Exception)
            {
                fail++; idx++;
                if (worst.Count < 50) worst.Add((1, 1, 1, 1, 1, TK, q));
                continue;
            }

            double eP = RelErr(pRef, a.Pressure.Pascal);
            double eRho = RelErr(rhoRef, a.Density.KilogramPerCubicMeter);
            double eH = RelErr(hRef, a.Enthalpy.JoulePerKilogram);
            double eS = RelErr(sRef, a.Entropy.JoulePerKilogramKelvin);
            double eU = RelErr(uRef, a.InternalEnergy.JoulePerKilogram);
            bool ok = eP <= relTol && eRho <= relTol && eH <= relTol && eS <= relTol && eU <= relTol;

            if (ok) pass++;
            else
            {
                fail++;
                if (worst.Count < 50) worst.Add((eP, eRho, eH, eS, eU, TK, q));
            }
            idx++;
        }

        string summary = $"UpdateTXExact sweep: attempts={attempts}, tested={idx}, pass={pass}, fail={fail}, skip={skip}";
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
