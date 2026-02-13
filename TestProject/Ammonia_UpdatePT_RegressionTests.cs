using Microsoft.VisualStudio.TestTools.UnitTesting;
using EngineeringUnits;
using EngineeringFluids.Fluids;
using static EngineeringFluids.Helmholtz.Phase;

[TestClass]
public class Ammonia_UpdatePT_RegressionTests
{
    //[TestMethod]
    //public void UpdatePT_StrictPhaseHint_Throws_WhenPhaseImpossible_But_AutoSucceeds()
    //{
    //    // Pick a definitely subcooled liquid point:
    //    // Use SharpFluids to get Tsat(P) and then go below it by a margin.
    //    // This ensures "Gas" is not a valid single-phase solution at that (P,T).
    //    var oracle = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);

    //    Pressure P = Pressure.FromBar(10);
    //    oracle.Pressure = P;

    //    var TsatObj = oracle.Tsat;
    //    if (oracle.FailState || TsatObj == null)
    //        Assert.Inconclusive("Oracle could not provide Tsat(P) for test setup.");

    //    double TsatK = TsatObj.Kelvin;

    //    // 30 K below saturation -> solidly subcooled liquid
    //    Temperature T = Temperature.FromKelvin(TsatK - 30.0);

    //    // 1) AUTO solver should succeed (it should pick liquid branch)
    //    var a_auto = new Ammonia();
    //    a_auto.UpdatePT(P, T);
    //    Assert.IsTrue(a_auto.Pressure.Pascal > 0, "Auto UpdatePT did not produce a valid state.");

    //    // Optional sanity: should be liquid (your Phase classifier may return Unknown if density not set,
    //    // but after UpdatePT it should have density)
    //    Assert.IsTrue(a_auto.Phase == Phases.Liquid || a_auto.Phase == Phases.Twophase,
    //        $"Auto UpdatePT unexpected phase: {a_auto.Phase}");

    //    // 2) STRICT gas solver must throw (no silent fallback allowed)
    //    var a_strict = new Ammonia();

    //    Assert.ThrowsException<InvalidOperationException>(
    //        () => a_strict.UpdatePT(P, T, Phases.Gas),
    //        "Strict UpdatePT(P,T,Gas) should throw for a subcooled liquid point, but it did not.");

    //    // 3) STRICT liquid solver should succeed
    //    var a_liq = new Ammonia();
    //    a_liq.UpdatePT(P, T, Phases.Liquid);
    //    Assert.IsTrue(a_liq.Pressure.Pascal > 0, "Strict UpdatePT(P,T,Liquid) did not produce a valid state.");

    //    // And it should match the AUTO solution reasonably (they should converge to same liquid root)
    //    double rho_auto = a_auto.Density!.KilogramPerCubicMeter;
    //    double rho_liq = a_liq.Density!.KilogramPerCubicMeter;

    //    double rel = System.Math.Abs(rho_auto - rho_liq) / System.Math.Max(System.Math.Abs(rho_auto), 1.0);
    //    Assert.IsTrue(rel < 1e-8, $"Auto vs strict-liquid density mismatch too large: rel={rel}");
    //}

    //[TestMethod]
    //public void UpdatePT_StrictPhaseHint_Throws_WhenPhaseImpossible_VaporCase()
    //{
    //    // Mirror case: pick a definitely superheated vapor point and ensure strict liquid throws.
    //    var oracle = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);

    //    Pressure P = Pressure.FromBar(10);
    //    oracle.Pressure = P;

    //    var TsatObj = oracle.Tsat;
    //    if (oracle.FailState || TsatObj == null)
    //        Assert.Inconclusive("Oracle could not provide Tsat(P) for test setup.");

    //    double TsatK = TsatObj.Kelvin;

    //    // 30 K above saturation -> solidly superheated vapor
    //    Temperature T = Temperature.FromKelvin(TsatK + 30.0);

    //    // AUTO solver should succeed
    //    var a_auto = new Ammonia();
    //    a_auto.UpdatePT(P, T);
    //    Assert.IsTrue(a_auto.Pressure.Pascal > 0, "Auto UpdatePT did not produce a valid state.");

    //    // STRICT liquid must throw
    //    var a_strict = new Ammonia();
    //    Assert.ThrowsException<InvalidOperationException>(
    //        () => a_strict.UpdatePT(P, T, Phases.Liquid),
    //        "Strict UpdatePT(P,T,Liquid) should throw for a superheated vapor point, but it did not.");

    //    // STRICT vapor should succeed
    //    var a_vap = new Ammonia();
    //    a_vap.UpdatePT(P, T, Phases.Gas);
    //    Assert.IsTrue(a_vap.Pressure.Pascal > 0, "Strict UpdatePT(P,T,Gas) did not produce a valid state.");
    //}
}