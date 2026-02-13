//using Microsoft.VisualStudio.TestTools.UnitTesting;
//using EngineeringUnits;
//using System;
//using System.Collections.Generic;
//using EngineeringFluids.Helmholtz.Solvers;

//[TestClass]
//public class CoolPropOracle_UpdateTX_Tests
//{
//    private static double RelErr(double expected, double actual)
//    {
//        double denom = Math.Max(Math.Abs(expected), 1.0);
//        return Math.Abs(actual - expected) / denom;
//    }

//    private static void SkipIfOracleFailed(SharpFluids.Fluid f, string context)
//    {
//        if (f.FailState)
//            Assert.Inconclusive($"SharpFluids oracle failed (FailState=true): {context}");
//    }

//    public static IEnumerable<object[]> TXPoints_Safe()
//    {
//        // Choose temperatures comfortably between triple and critical for NH3
//        // (You can expand later once stable.)
//        double[] tempsK = { 210, 220, 240, 260, 280, 300, 330, 360, 380, 395 };
//        double[] qs = { 0.05, 0.2, 0.5, 0.8, 0.95 };

//        foreach (var T in tempsK)
//            foreach (var q in qs)
//                yield return new object[] { T, q };
//    }

//    [DataTestMethod]
//    [DynamicData(nameof(TXPoints_Safe), DynamicDataSourceType.Method)]
//    public void UpdateTX_Matches_SharpFluids(double tK, double q)
//    {
//        const double absTolT_K = 5e-3;     // T is input → should be tight
//        const double relTolH = 1e-3;       // 0.1%
//        const double relTolRho = 5e-3;     // 0.5%
//        const double absTolQ = 1e-12;

//        // Loose oracle pressure tolerance (correlation vs EOS)
//        const double absTolP_oracle_Pa = 50_000; // 0.5 bar

//        Temperature T = Temperature.FromKelvin(tK);

//        // Oracle (CoolProp via SharpFluids)
//        var refFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        refFluid.UpdateXT(q, T);
//        SkipIfOracleFailed(refFluid, $"UpdateXT oracle at T={tK}K q={q}");

//        double p_ref = refFluid.Pressure!.Pascal;
//        double rho_ref = refFluid.Density!.KilogramPerCubicMeter;
//        double h_ref = refFluid.Enthalpy!.JoulePerKilogram;
//        double T_ref = refFluid.Temperature!.Kelvin;
//        double q_ref = refFluid.Quality;

//        // Our implementation
//        var a = new EngineeringFluids.Fluids.Ammonia();
//        a.UpdateTX(T, q);

//        // Contract: T and q exact
//        Assert.IsTrue(Math.Abs(a.Temperature!.Kelvin - T_ref) < absTolT_K, "T mismatch vs oracle");
//        Assert.IsTrue(Math.Abs(a.Quality - q) < absTolQ, "Quality mismatch (contract)");
//        Assert.IsTrue(Math.Abs(q_ref - q) < 1e-9, "Oracle quality mismatch (unexpected)");

//        // Contract: Psat must match OUR correlation tightly
//        var psat_expected = EngineeringFluids.Helmholtz.Saturation.CalculateSaturationPressure(T);
//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - psat_expected.Pascal) < 1.0,
//            $"Psat mismatch vs own correlation: expected={psat_expected.Pascal}, actual={a.Pressure.Pascal}");

//        // Oracle sanity: Psat roughly close to CoolProp (loose)
//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - p_ref) < absTolP_oracle_Pa,
//            $"Psat mismatch vs oracle (loose): oracle={p_ref}, actual={a.Pressure.Pascal}");

//        // Properties: modest tolerances
//        Assert.IsTrue(RelErr(h_ref, a.Enthalpy.SI) < relTolH, "h mismatch vs oracle");
//        Assert.IsTrue(RelErr(rho_ref, a.Density!.KilogramPerCubicMeter) < relTolRho, "rho mismatch vs oracle");
//    }

//    [TestMethod]
//    public void UpdateTX_RoundTrip_To_UpdatePX_SelfConsistency()
//    {
//        // TX uses correlation saturation; PX uses EOS SolveAtP (usually).
//        // Expect very small differences in endpoints => tiny enthalpy deltas.
//        const double relTol = 1e-7; // was 1e-10 (too strict)

//        Temperature T = Temperature.FromKelvin(280);
//        double q = 0.4;

//        var a1 = new EngineeringFluids.Fluids.Ammonia();
//        a1.UpdateTX(T, q);
//        double p = a1.Pressure.Pascal;
//        double h1 = a1.Enthalpy.SI;

//        var a2 = new EngineeringFluids.Fluids.Ammonia();
//        a2.UpdatePX(Pressure.FromPascal(p), q);
//        double h2 = a2.Enthalpy.SI;

//        double rel = Math.Abs(h1 - h2) / Math.Max(Math.Abs(h1), 1.0);
//        Assert.IsTrue(rel < relTol, $"TX->PX self-consistency failed: h1={h1}, h2={h2}, rel={rel}");
//    }
//}