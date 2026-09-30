using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;
using System.Threading.Tasks;
using Phases = EngineeringFluids.Helmholtz.Phase.Phases;

namespace TestProject;

// Bugs found while matching SharpFluids, and the behaviour that replaced them.
[TestClass]
public class RegressionTests
{
    // UpdatePH just above Pc silently returned a wrong state: a warm-started liquid density solve converged on a
    // spurious root of the EOS's unstable loop, and the bisection then accepted a collapsed bracket without
    // checking the residual (h = 1.12 MJ/kg came back as 0.865 MJ/kg at 375.5 K).
    [TestMethod]
    public void UpdatePH_JustAboveTheCriticalPressure_FindsTheRequestedEnthalpy()
    {
        foreach (double h in new[] { 1.119620e6, 1.2e6, 1.203165e6 })
        {
            var a = new Ammonia();
            a.UpdatePH(Pressure.FromBar(121.55), Enthalpy.FromJoulePerKilogram(h));
            Assert.AreEqual(h, a.Enthalpy.JoulePerKilogram, 1e-6 * h, $"h={h}");
            Assert.IsTrue(a.Temperature.Kelvin > 405, $"T={a.Temperature.Kelvin} K is the old wrong liquid state");
        }
    }

    // Every state PH/PS return must match the input: they either find it or throw, never return something else.
    [TestMethod]
    public void UpdatePHAndPS_NeverReturnAStateThatMissesTheInput()
    {
        foreach (double p in Enumerable.Range(0, 30).Select(i => 1e5 * Math.Pow(1000.0, i / 29.0)))
        {
            foreach (double h in Enumerable.Range(0, 40).Select(i => -5e4 + i * 2.2e6 / 39))
            {
                var a = new Ammonia();
                try { a.UpdatePH(Pressure.FromPascal(p), Enthalpy.FromJoulePerKilogram(h)); }
                catch (InvalidOperationException) { continue; }

                // Within 50 J/kg of the dome UpdatePH snaps to the saturation line by design (UpdatePHExact does not)
                Assert.AreEqual(h, a.Enthalpy.JoulePerKilogram, a.Quality is 0 or 1 ? 50.0 : 1e-6 * Math.Max(Math.Abs(h), 1), $"PH p={p} Pa h={h}");
            }
            foreach (double s in Enumerable.Range(0, 40).Select(i => -200 + i * 7800.0 / 39))
            {
                var a = new Ammonia();
                try { a.UpdatePS(Pressure.FromPascal(p), SpecificEntropy.FromJoulePerKilogramKelvin(s)); }
                catch (InvalidOperationException) { continue; }

                Assert.AreEqual(s, a.Entropy.JoulePerKilogramKelvin, a.Quality is 0 or 1 ? 1e-3 : 1e-6 * Math.Max(Math.Abs(s), 1), $"PS p={p} Pa s={s}");
            }
        }
    }

    // SharpFluids falls back to exactly LimitTemperatureMin (= the triple point) in AddPower/RemovePower; UpdatePT
    // used to fail there at every pressure because the liquid density guess skipped the triple point itself.
    [TestMethod]
    public void UpdatePT_AtTheTripleTemperature_Works()
    {
        foreach (double p in new[] { 1e4, 1e5, 1e6, 1e7, 1e8 })
        {
            var a = new Ammonia();
            a.UpdatePT(Pressure.FromPascal(p), a.LimitTemperatureMin);
            Assert.AreEqual(p, a.Pressure.Pascal, 1e-6 * p);
            Assert.IsTrue(a.Density.KilogramPerCubicMeter > 700, "liquid");
        }
    }

    // In the last ~0.3 K below Tc the saturation Newton, started from the old float ancillaries, slid onto rhoL == rhoV
    // (which also satisfies equal pressure and fugacity) or onto a nearby spurious root: from 405.24 K to 405.38 K it
    // returned rhoL/rhoV = 1.0002 where the dome is 1.378 wide, p 1.2e-4 off, and it threw from 405.4 K. The Gao dome
    // continues to Tc, and SolveAtT must now follow it all the way.
    [TestMethod]
    public void SolveAtT_FollowsTheRealDomeUpToTc()
    {
        var a = new Ammonia();
        foreach (var d in AncillaryTests.DomeReference)
        {
            var sat = a.SolveAtT(Temperature.FromKelvin(d.T));
            Assert.AreEqual(d.P, sat.Psat.Pascal, 1e-9 * d.P, $"p at {d.T} K");
            Assert.AreEqual(d.RhoL, sat.RhomolarL.MolesPerCubicMeter, 1e-7 * d.RhoL, $"rhoL at {d.T} K");
            Assert.AreEqual(d.RhoV, sat.RhomolarV.MolesPerCubicMeter, 1e-7 * d.RhoV, $"rhoV at {d.T} K");

            var byP = a.SolveAtP(Pressure.FromPascal(d.P));
            Assert.AreEqual(d.T, byP.T.Kelvin, 1e-8, $"SolveAtP at {d.T} K");
        }

        // and it never returns a trivial or spurious root, down to 10 microK below Tc
        foreach (double t in Enumerable.Range(0, 60).Select(i => 400.0 + i * (405.55999 - 400.0) / 59))
        {
            var sat = a.SolveAtT(Temperature.FromKelvin(t));
            double w = Math.Log((double)(sat.RhomolarL / sat.RhomolarV));
            double wAncillary = Math.Log((double)(BubbleDensity.Density(Temperature.FromKelvin(t)) / DewDensity.Density(Temperature.FromKelvin(t))));
            Assert.AreEqual(wAncillary, w, 1e-3 * wAncillary, $"dome width at {t} K");
        }
    }

    // The exact two-phase flashes solve the EOS, so they agree with each other to the solver tolerance
    [TestMethod]
    public void ExactTwoPhaseFlashes_AreConsistentWithEachOther()
    {
        foreach (double t in Enumerable.Range(0, 20).Select(i => 200.0 + i * (404.0 - 200.0) / 19))
        {
            var byT = new Ammonia();
            byT.UpdateTXExact(Temperature.FromKelvin(t), 0.4);

            var byP = new Ammonia();
            byP.UpdatePXExact(byT.Pressure, 0.4);
            Assert.AreEqual(t, byP.Temperature.Kelvin, 1e-8 * t, $"T at {t} K");

            var byH = new Ammonia();
            byH.UpdatePHExact(byT.Pressure, byT.Enthalpy);
            Assert.AreEqual(0.4, byH.Quality, 1e-7, $"q from PH at {t} K");

            var byD = new Ammonia();
            byD.UpdateDT(byT.Density, byT.Temperature);
            Assert.AreEqual(0.4, byD.Quality, 1e-7, $"q from DT at {t} K");
        }
    }

    [TestMethod]
    public void UpdatePXExact_AtTheCriticalPressure_GivesTheCriticalPoint()
    {
        var a = new Ammonia();
        a.UpdatePXExact(a.CriticalPressure, 1.0);
        Assert.AreEqual(Phases.CriticalPoint, a.Phase);
        Assert.AreEqual(a.CriticalTemperature.Kelvin, a.Temperature.Kelvin);
        Assert.AreEqual(233.25, a.Density.KilogramPerCubicMeter, 1e-2);
        Assert.AreEqual(1.0, a.Quality);
    }

    // Phase labels follow CoolProp: above Pc by T vs Tc, above Tc and below Pc SupercriticalGas
    [TestMethod]
    public void Phase_FollowsCoolPropsLabels()
    {
        (double p, double t, Phases expected)[] cases =
        {
            (50e5, 450, Phases.SupercriticalGas),
            (150e5, 350, Phases.SupercriticalLiquid),
            (150e5, 450, Phases.Supercritical),
            (5e5, 300, Phases.Gas),
            (50e5, 250, Phases.Liquid),
        };
        foreach (var (p, t, expected) in cases)
        {
            var a = new Ammonia();
            a.UpdatePT(Pressure.FromPascal(p), Temperature.FromKelvin(t));
            Assert.AreEqual(expected, a.Phase, $"PT {p} Pa, {t} K");

            var b = new Ammonia();
            b.UpdatePHExact(Pressure.FromPascal(p), a.Enthalpy);
            Assert.AreEqual(expected, b.Phase, $"PH {p} Pa, {t} K");
        }
    }

    [TestMethod]
    public void SaturationHelpers_BehaveLikeSharpFluids()
    {
        var two = new Ammonia();
        two.UpdatePXExact(Pressure.FromBar(10), 0.3);
        Assert.AreEqual(two.Temperature.Kelvin, two.Tsat.Kelvin);
        Assert.IsTrue(two.LiquidDensity > two.Density && two.Density > two.GasDensity);

        var gas = new Ammonia();
        gas.UpdatePT(Pressure.FromBar(10), Temperature.FromKelvin(350));
        Assert.AreEqual(two.Temperature.Kelvin, gas.Tsat.Kelvin, 1e-6, "Tsat at the same pressure");
        Assert.AreEqual(gas.Density, gas.GasDensity);
        Assert.AreEqual(gas.Density, gas.LiquidDensity);

        var above = new Ammonia();
        above.UpdatePT(Pressure.FromBar(200), Temperature.FromKelvin(450));
        Assert.AreEqual(above.CriticalTemperature.Kelvin, above.Tsat.Kelvin);

        var below = new Ammonia();
        below.UpdatePT(Pressure.FromPascal(3000), Temperature.FromKelvin(250));
        Assert.IsTrue(double.IsNaN(below.Tsat.Kelvin), "no saturation state below the triple-point pressure");
    }

    // Once a state is set, many threads may read the same instance: the cache is published safely
    [TestMethod]
    public void ConcurrentReads_OfOneInstance_AllSeeTheSameValues()
    {
        for (int round = 0; round < 50; round++)
        {
            var shared = new Ammonia();
            shared.UpdatePXExact(Pressure.FromBar(5 + round), 0.4); // fresh, unread cache every round

            var reference = new Ammonia();
            reference.UpdatePXExact(Pressure.FromBar(5 + round), 0.4);
            double[] expected = Read(reference);

            var results = new double[16][];
            Parallel.For(0, results.Length, i => results[i] = Read(shared));

            foreach (double[] r in results)
                CollectionAssert.AreEqual(expected, r, $"round {round}");
        }

        static double[] Read(Ammonia a) => new[]
        {
            a.Enthalpy.SI, a.Entropy.SI, a.Cp.SI, a.Cv.SI, a.DynamicViscosity.SI, a.Conductivity.SI, a.Prandtl,
            a.SurfaceTension.SI, a.Tsat.SI, a.GasDensity.SI, a.LiquidDensity.SI, a.AlphaR_dDelta, a.Compressibility,
        };
    }
}
