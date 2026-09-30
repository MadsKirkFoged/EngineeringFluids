using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace TestProject;

// Ammonia caches the Helmholtz derivatives and transport properties per (Temperature, Density) state.
// Every way the state can change must invalidate that cache: after any change, every property has to be
// bit-identical to a fresh instance put directly into the new state.
[TestClass]
public class StateCacheTests
{
    private static void AssertSameState(Ammonia expected, Ammonia actual)
    {
        Assert.AreEqual(expected.Temperature.Kelvin, actual.Temperature.Kelvin);
        Assert.AreEqual(expected.Density.KilogramPerCubicMeter, actual.Density.KilogramPerCubicMeter);
        Assert.AreEqual(expected.Pressure.Pascal, actual.Pressure.Pascal);
        Assert.AreEqual(expected.Enthalpy.JoulePerKilogram, actual.Enthalpy.JoulePerKilogram);
        Assert.AreEqual(expected.Entropy.JoulePerKilogramKelvin, actual.Entropy.JoulePerKilogramKelvin);
        Assert.AreEqual(expected.InternalEnergy.JoulePerKilogram, actual.InternalEnergy.JoulePerKilogram);
        Assert.AreEqual(expected.Cp.JoulePerKilogramKelvin, actual.Cp.JoulePerKilogramKelvin);
        Assert.AreEqual(expected.Cv.JoulePerKilogramKelvin, actual.Cv.JoulePerKilogramKelvin);
        Assert.AreEqual(expected.SoundSpeed.MeterPerSecond, actual.SoundSpeed.MeterPerSecond);
        Assert.AreEqual(expected.DynamicViscosity.PascalSecond, actual.DynamicViscosity.PascalSecond);
        Assert.AreEqual(expected.Conductivity.WattPerMeterKelvin, actual.Conductivity.WattPerMeterKelvin);
        Assert.AreEqual(expected.SurfaceTension.NewtonPerMeter, actual.SurfaceTension.NewtonPerMeter);
        Assert.AreEqual(expected.Quality, actual.Quality);
        Assert.AreEqual(expected.Phase, actual.Phase);
        Assert.AreEqual(expected.AlphaR_dDelta, actual.AlphaR_dDelta);
    }

    // Reads every cached property so the cache is warm before the state changes
    private static void ReadEverything(Ammonia a) => AssertSameState(a, a);

    [TestMethod]
    public void SettingTemperature_RefreshesEveryProperty()
    {
        var a = new Ammonia { Temperature = Temperature.FromKelvin(400), Density = Density.FromKilogramPerCubicMeter(10) };
        ReadEverything(a);

        a.Temperature = Temperature.FromKelvin(450);

        AssertSameState(new Ammonia { Temperature = Temperature.FromKelvin(450), Density = Density.FromKilogramPerCubicMeter(10) }, a);
    }

    [TestMethod]
    public void SettingDensity_RefreshesEveryProperty()
    {
        var a = new Ammonia { Temperature = Temperature.FromKelvin(400), Density = Density.FromKilogramPerCubicMeter(10) };
        ReadEverything(a);

        a.Density = Density.FromKilogramPerCubicMeter(20);

        AssertSameState(new Ammonia { Temperature = Temperature.FromKelvin(400), Density = Density.FromKilogramPerCubicMeter(20) }, a);
    }

    [TestMethod]
    public void SinglePhaseThenTwoPhase_MatchesAFreshTwoPhaseState()
    {
        var a = new Ammonia();
        a.UpdatePT(Pressure.FromBar(20), Temperature.FromKelvin(400));
        ReadEverything(a);

        a.UpdatePX(Pressure.FromBar(10), 0.3);

        var fresh = new Ammonia();
        fresh.UpdatePX(Pressure.FromBar(10), 0.3);
        AssertSameState(fresh, a);
    }

    [TestMethod]
    public void TwoPhaseThenSinglePhase_MatchesAFreshSinglePhaseState()
    {
        var a = new Ammonia();
        a.UpdatePX(Pressure.FromBar(10), 0.3);
        ReadEverything(a);

        a.UpdatePT(Pressure.FromBar(20), Temperature.FromKelvin(400));

        var fresh = new Ammonia();
        fresh.UpdatePT(Pressure.FromBar(20), Temperature.FromKelvin(400));
        AssertSameState(fresh, a);
    }

    [TestMethod]
    public void TwoPhaseThenAnotherTwoPhase_DropsTheOldSaturatedStates()
    {
        var a = new Ammonia();
        a.UpdatePX(Pressure.FromBar(10), 0.0);
        ReadEverything(a);

        a.UpdateTX(Temperature.FromKelvin(250), 1.0);

        var fresh = new Ammonia();
        fresh.UpdateTX(Temperature.FromKelvin(250), 1.0);
        AssertSameState(fresh, a);
    }
}
