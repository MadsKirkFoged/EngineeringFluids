using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using static EngineeringFluids.Helmholtz.Phase;

namespace TestProject;

// Cp, Cv, speed of sound, viscosity, thermal conductivity, Prandtl and surface tension against CoolProp.
//
// CoolProp is set to the SAME (T, rho) state EngineeringFluids computed, so these tests check the property
// formulas and correlations themselves, independent of how accurately the flash found that state.
//
// The oracle is CoolProp's own PropsSI, called straight into the CoolProp64 library that the SharpFluids package
// ships - not SharpFluids' AbstractState wrapper, which returns 0 for the speed of sound whenever 0 <= Q <= 1 and
// hides the surface tension at Q = 0 and Q = 1.
[TestClass]
public class PropertyTests
{
    private static class CoolProp
    {
        [DllImport("CoolProp64", EntryPoint = "CSharp_PropsSI")]
        public static extern double PropsSI(string output, string name1, double value1, string name2, double value2, string fluid);

        public static double AtTD(string output, Ammonia a)
            => PropsSI(output, "T", a.Temperature.Kelvin, "Dmass", a.Density.KilogramPerCubicMeter, "Ammonia");
    }

    private static double RelErr(double actual, double expected) => Math.Abs(actual - expected) / Math.Abs(expected);

    private static IEnumerable<Ammonia> SinglePhaseGrid()
    {
        foreach (double p in Enumerable.Range(0, 20).Select(i => Math.Exp(Math.Log(1e4) + i * (Math.Log(1e8) - Math.Log(1e4)) / 19)))
        {
            foreach (double t in Enumerable.Range(0, 20).Select(i => 196.0 + i * (725.0 - 196.0) / 19))
            {
                var a = new Ammonia();
                try
                {
                    a.UpdatePT(Pressure.FromPascal(p), Temperature.FromKelvin(t));
                }
                catch (InvalidOperationException)
                {
                    continue; // (T,P) on the saturation line: UpdatePT refuses those by design
                }
                yield return a;
            }
        }
    }

    [TestMethod]
    public void SinglePhase_CaloricAndTransportProperties_MatchCoolProp()
    {
        const double relTol = 1e-10;
        var worst = new Dictionary<string, (double err, string at)>();

        void Check(string name, double actual, double expected, Ammonia a)
        {
            double err = RelErr(actual, expected);
            if (!worst.TryGetValue(name, out var w) || !(err <= w.err))
                worst[name] = (err, $"T={a.Temperature.Kelvin:G6} K, rho={a.Density.KilogramPerCubicMeter:G6} kg/m3");
        }

        int count = 0;
        foreach (Ammonia a in SinglePhaseGrid())
        {
            count++;
            Check("Cp", a.Cp.JoulePerKilogramKelvin, CoolProp.AtTD("Cpmass", a), a);
            Check("Cv", a.Cv.JoulePerKilogramKelvin, CoolProp.AtTD("Cvmass", a), a);
            Check("SoundSpeed", a.SoundSpeed.MeterPerSecond, CoolProp.AtTD("A", a), a);
            Check("DynamicViscosity", a.DynamicViscosity.PascalSecond, CoolProp.AtTD("V", a), a);
            Check("Conductivity", a.Conductivity.WattPerMeterKelvin, CoolProp.AtTD("L", a), a);
            Check("Prandtl", a.Prandtl, CoolProp.AtTD("Prandtl", a), a);
        }

        Assert.IsTrue(count > 350, $"Only {count} grid states - the grid no longer covers the domain.");
        foreach (var (name, (err, at)) in worst)
            Assert.IsTrue(err <= relTol, $"{name}: relative error {err:E2} at {at}");
    }

    // Matches the CoolProp values SharpFluids' own AmmoniaTests.UpdateDH asserts for rho = 50 kg/m3,
    // h = 2148349.72016398 J/kg (supercritical, 278.39 C, 117.97 bar).
    [TestMethod]
    public void SupercriticalPoint_MatchesSharpFluidsReferenceValues()
    {
        var a = new Ammonia();
        a.UpdatePT(Pressure.FromBar(117.968517492664), Temperature.FromDegreeCelsius(278.392108051419));

        Assert.AreEqual(50, a.Density.KilogramPerCubicMeter, 1e-6);
        Assert.AreEqual(3251.34163624504, a.Cp.JoulePerKilogramKelvin, 1e-5);
        Assert.AreEqual(2270.36197650328, a.Cv.JoulePerKilogramKelvin, 1e-5);
        Assert.AreEqual(544.366706172041, a.SoundSpeed.MeterPerSecond, 1e-6);
        Assert.AreEqual(2.0869708899873E-05, a.DynamicViscosity.PascalSecond, 1e-13);
        Assert.AreEqual(0.0722711120542153, a.Conductivity.WattPerMeterKelvin, 1e-10);
        Assert.AreEqual(0.93888901877652886, a.Prandtl, 1e-8);
        Assert.IsTrue(double.IsNaN(a.SurfaceTension.NewtonPerMeter), "Surface tension is only defined in two-phase.");
    }

    // In two-phase CoolProp evaluates cp, cv and the transport properties at the bulk state (saturation T, mixture
    // density). That point is inside the dome and has no physical meaning, but it is what CoolProp - and so every
    // SharpFluids user today - gets, so the same bulk convention is kept.
    [TestMethod]
    public void TwoPhase_BulkProperties_MatchCoolPropAtTheSameMixtureDensity()
    {
        const double relTol = 1e-10;
        foreach (double t in Enumerable.Range(0, 20).Select(i => 200.0 + i * (400.0 - 200.0) / 19))
        {
            foreach (double q in new[] { 0.1, 0.5, 0.9 })
            {
                var a = new Ammonia();
                a.UpdateTX(Temperature.FromKelvin(t), q);
                string at = $"T={t:G6} K, q={q}";

                Assert.IsTrue(RelErr(a.Cp.JoulePerKilogramKelvin, CoolProp.AtTD("Cpmass", a)) <= relTol, $"Cp at {at}");
                Assert.IsTrue(RelErr(a.Cv.JoulePerKilogramKelvin, CoolProp.AtTD("Cvmass", a)) <= relTol, $"Cv at {at}");
                Assert.IsTrue(RelErr(a.DynamicViscosity.PascalSecond, CoolProp.AtTD("V", a)) <= relTol, $"viscosity at {at}");
                Assert.IsTrue(RelErr(a.Conductivity.WattPerMeterKelvin, CoolProp.AtTD("L", a)) <= relTol, $"conductivity at {at}");
            }
        }
    }

    [TestMethod]
    public void SurfaceTension_MatchesCoolPropInTwoPhase()
    {
        foreach (double t in Enumerable.Range(0, 40).Select(i => 196.0 + i * (404.5 - 196.0) / 39))
        {
            foreach (double q in new[] { 0.0, 0.5, 1.0 })
            {
                var a = new Ammonia();
                a.UpdateTX(Temperature.FromKelvin(t), q);
                double expected = CoolProp.PropsSI("I", "T", t, "Q", q, "Ammonia");
                Assert.AreEqual(expected, a.SurfaceTension.NewtonPerMeter, 1e-15 * Math.Abs(expected), $"T={t:G6} K, q={q}");
            }
        }
    }

    // Like CoolProp: on the saturation lines it is the saturated liquid/vapor value, in between it is undefined
    // (it depends on how the phases are distributed). The saturated densities come from the ancillary fits, whose
    // error (largest within ~10 K of Tc) is what limits the agreement here, not the speed-of-sound formula.
    [TestMethod]
    public void SoundSpeed_InTwoPhase_IsOnlyDefinedOnTheSaturationLines()
    {
        foreach (double t in Enumerable.Range(0, 40).Select(i => 196.0 + i * (404.5 - 196.0) / 39))
        {
            var mid = new Ammonia();
            mid.UpdateTX(Temperature.FromKelvin(t), 0.5);
            Assert.IsTrue(double.IsNaN(mid.SoundSpeed.MeterPerSecond), $"q=0.5 at T={t:G6} K should be NaN");

            foreach (double q in new[] { 0.0, 1.0 })
            {
                var a = new Ammonia();
                a.UpdateTX(Temperature.FromKelvin(t), q);
                double expected = CoolProp.PropsSI("A", "T", t, "Q", q, "Ammonia");
                Assert.IsTrue(RelErr(a.SoundSpeed.MeterPerSecond, expected) <= 2e-4, $"T={t:G6} K, q={q}: {a.SoundSpeed.MeterPerSecond} vs {expected}");
            }
        }
    }

    [TestMethod]
    public void SinglePhase_Phase_IsGasOrLiquid_AndSurfaceTensionIsNaN()
    {
        var a = new Ammonia();
        a.UpdatePT(Pressure.FromBar(5), Temperature.FromKelvin(300));
        Assert.AreEqual(Phases.Gas, a.Phase);
        Assert.IsTrue(double.IsNaN(a.SurfaceTension.NewtonPerMeter));
    }
}
