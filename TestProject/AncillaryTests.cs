using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Linq;
using System.Runtime.InteropServices;

namespace TestProject;

// The saturation ancillary tables (SaturationPressure, SaturationTemperature, BubbleDensity, DewDensity) against CoolProp,
// and above 405 K - where CoolProp's saturation solver breaks down - against the Gao EOS dome itself.
[TestClass]
public class AncillaryTests
{
    private static class CoolProp
    {
        [DllImport("CoolProp64", EntryPoint = "CSharp_PropsSI")]
        public static extern double PropsSI(string output, string name1, double value1, string name2, double value2, string fluid);
    }

    private const double Tc = 405.56;
    private const double Ttriple = 195.495;
    private const double PcEos = 11363391.162815318; // the EOS's own pressure at (Tc, rhoc)

    // (T K, p Pa, rhoL, rhoV mol/m3) on the EOS dome, followed by continuation (2x2 Newton on equal p and g) from 405 K.
    // 405.24 K is where SolveAtT used to return a spurious near-trivial root.
    internal static readonly (double T, double P, double RhoL, double RhoV)[] DomeReference =
    {
        (405.24, 11299006.227379337, 15887.841961359507, 11530.659149434367),
        (405.3, 11311044.066892788, 15728.91767891877, 11684.590774476324),
        (405.4, 11331141.313654043, 15407.906088896994, 11997.925081678435),
        (405.5, 11351282.853309637, 14928.144762649683, 12471.796999831826),
        (405.55, 11361371.636245389, 14405.652986763143, 12993.50521998949),
        (405.559, 11363189.169488676, 14051.60768873159, 13347.362902230341),
    };

    private static double Rel(double expected, double actual) => Math.Abs(actual / expected - 1);

    [TestMethod]
    public void Tables_MatchCoolProp_FromTheTriplePointTo405K()
    {
        double worstP = 0, worstL = 0, worstV = 0, worstT = 0;
        foreach (double t in Enumerable.Range(0, 400).Select(i => Ttriple + i * (405.0 - Ttriple) / 399))
        {
            double p = CoolProp.PropsSI("P", "T", t, "Q", 0, "Ammonia");
            double rhoL = CoolProp.PropsSI("Dmolar", "T", t, "Q", 0, "Ammonia");
            double rhoV = CoolProp.PropsSI("Dmolar", "T", t, "Q", 1, "Ammonia");
            var T = Temperature.FromKelvin(t);

            worstP = Math.Max(worstP, Rel(p, SaturationPressure.Pressure(T).Pascal));
            worstL = Math.Max(worstL, Rel(rhoL, BubbleDensity.Density(T).MolesPerCubicMeter));
            worstV = Math.Max(worstV, Rel(rhoV, DewDensity.Density(T).MolesPerCubicMeter));
            worstT = Math.Max(worstT, Math.Abs(SaturationTemperature.Temperature(Pressure.FromPascal(p)).Kelvin - t));
        }
        // measured over 51,000 points when the tables were fitted: 1.0e-10, 2.0e-9, 9.8e-11 and 1.2e-8 K
        Assert.IsTrue(worstP < 3e-10, $"p: {worstP:E2}");
        Assert.IsTrue(worstL < 5e-9, $"rhoL: {worstL:E2}");
        Assert.IsTrue(worstV < 3e-10, $"rhoV: {worstV:E2}");
        Assert.IsTrue(worstT < 3e-8, $"T(p): {worstT:E2} K");
    }

    [TestMethod]
    public void Tables_FollowTheEosDome_AboveCoolPropsRange()
    {
        foreach (var d in DomeReference)
        {
            var T = Temperature.FromKelvin(d.T);
            Assert.IsTrue(Rel(d.P, SaturationPressure.Pressure(T).Pascal) < 1e-9, $"p at {d.T} K");
            Assert.IsTrue(Rel(d.RhoL, BubbleDensity.Density(T).MolesPerCubicMeter) < 1e-7, $"rhoL at {d.T} K");
            Assert.IsTrue(Rel(d.RhoV, DewDensity.Density(T).MolesPerCubicMeter) < 1e-7, $"rhoV at {d.T} K");
            Assert.AreEqual(d.T, SaturationTemperature.Temperature(Pressure.FromPascal(d.P)).Kelvin, 1e-8, $"T(p) at {d.T} K");
        }
    }

    // Bisection and Newton on these curves (UpdatePX, SolveAtP, the dome-temperature search in UpdateDP) need them
    // monotonic, and T(p) must invert p(T)
    [TestMethod]
    public void Tables_AreMonotonic_AndInvertEachOther()
    {
        double prevP = 0, prevL = double.MaxValue, prevV = 0;
        for (double t = Ttriple; t < Tc; t += 0.00037)
        {
            var T = Temperature.FromKelvin(t);
            double p = SaturationPressure.Pressure(T).Pascal, rhoL = BubbleDensity.Density(T).MolesPerCubicMeter, rhoV = DewDensity.Density(T).MolesPerCubicMeter;
            Assert.IsTrue(p > prevP && rhoL < prevL && rhoV > prevV && rhoL > rhoV, $"not monotonic at {t} K");
            Assert.AreEqual(t, SaturationTemperature.Temperature(Pressure.FromPascal(p)).Kelvin, 3e-8, $"T(p(T)) at {t} K");
            prevP = p; prevL = rhoL; prevV = rhoV;
        }
    }

    [TestMethod]
    public void Tables_EndAtTheCriticalPoint_AndClampOutsideTheirRange()
    {
        foreach (double t in new[] { Tc, Tc + 50 })
        {
            var T = Temperature.FromKelvin(t);
            Assert.AreEqual(PcEos, SaturationPressure.Pressure(T).Pascal, 1e-9 * PcEos);
            Assert.AreEqual(13696, BubbleDensity.Density(T).MolesPerCubicMeter, 1.0);
            Assert.AreEqual(13696, DewDensity.Density(T).MolesPerCubicMeter, 1.0);
        }
        Assert.AreEqual(Tc, SaturationTemperature.Temperature(Pressure.FromPascal(PcEos)).Kelvin, 1e-8);
        Assert.AreEqual(Tc, SaturationTemperature.Temperature(Pressure.FromBar(200)).Kelvin, 1e-8);

        // below the triple point: the triple-point values (its EOS saturation pressure is 6055.8 Pa)
        Assert.AreEqual(6055.81357399642, SaturationPressure.Pressure(Temperature.FromKelvin(150)).Pascal, 1e-5);
        Assert.AreEqual(Ttriple, SaturationTemperature.Temperature(Pressure.FromPascal(100)).Kelvin, 1e-8);

        Assert.IsTrue(double.IsNaN(SaturationPressure.Pressure(Temperature.FromKelvin(double.NaN)).Pascal));
        Assert.IsTrue(double.IsNaN(BubbleDensity.Density(Temperature.FromKelvin(double.NaN)).MolesPerCubicMeter));
        Assert.IsTrue(double.IsNaN(SaturationTemperature.Temperature(Pressure.FromPascal(double.NaN)).Kelvin));
    }
}
