using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace TestProject;

[TestClass]
public class DerivativeTests
{
    private static readonly Ammonia Ref = new();
    private static readonly Density RhoCritMass = Ref.Critical.MolarDensity * Ref.MolarMass;

    private static double AdaptiveEps(double delta)
    {
        double eps = 1e-6 * Math.Max(1.0, Math.Abs(delta));
        if (delta - eps <= 0.0)
            eps = 0.5 * delta; // keep delta-eps > 0
        return eps;
    }

    private static Ammonia StateAt(double T, double delta)
        => new Ammonia { Temperature = Temperature.FromKelvin(T), Density = RhoCritMass * delta };

    private static void AssertDDelta2MatchesFiniteDifference(double T, double delta, string label)
    {
        double eps = AdaptiveEps(delta);

        double num = (StateAt(T, delta + eps).AlphaR_dDelta - StateAt(T, delta - eps).AlphaR_dDelta) / (2.0 * eps);
        double ana = StateAt(T, delta).AlphaR_dDelta2;

        // If this fails badly, it usually means a sign error in one term.
        Assert.AreEqual(num, ana, 1e-7, $"{label} mismatch: num={num}, ana={ana}, delta={delta}, eps={eps}");
    }

    [TestMethod]
    public void AlphaR_dDelta2_MatchesFiniteDifference_VaporLikePoint()
    {
        // T = 400 K, rho_mass = 9 kg/m3 => delta ~ 0.0386
        AssertDDelta2MatchesFiniteDifference(400.0, 9.0 / RhoCritMass.KilogramPerCubicMeter, "Vapor-like point");
    }

    [TestMethod]
    public void AlphaR_dDelta2_MatchesFiniteDifference_LiquidLikePoint()
    {
        AssertDDelta2MatchesFiniteDifference(400.0, 1.4748450193131541, "Liquid-like point");
    }

    // ------------------------------------------------------------------
    // Second tau derivatives (used by Cp, Cv and the speed of sound): central differences of the
    // first derivatives in tau at fixed delta. Each point is (T, delta) - vapor, liquid, near-critical
    // and supercritical - so every residual class (power, Gaussian, Gao-B) contributes noticeably.
    // ------------------------------------------------------------------
    private static readonly double Tc = Ref.Critical.Temperature.Kelvin;

    private static Ammonia StateAtTau(double tau, double delta) => StateAt(Tc / tau, delta);

    public static IEnumerable<object[]> TauDerivativePoints()
    {
        yield return new object[] { 400.0, 0.0386 };  // vapor
        yield return new object[] { 280.0, 2.73 };    // liquid
        yield return new object[] { 406.0, 1.0 };     // near-critical
        yield return new object[] { 600.0, 0.5 };     // supercritical
    }

    [DataTestMethod]
    [DynamicData(nameof(TauDerivativePoints), DynamicDataSourceType.Method)]
    public void SecondTauDerivatives_MatchFiniteDifference(double T, double delta)
    {
        double tau = Tc / T;
        double h = 1e-6 * tau;
        Ammonia plus = StateAtTau(tau + h, delta), minus = StateAtTau(tau - h, delta), at = StateAtTau(tau, delta);

        double num0 = (plus.Alpha0_dTau - minus.Alpha0_dTau) / (2.0 * h);
        double numR = (plus.AlphaR_dTau - minus.AlphaR_dTau) / (2.0 * h);

        Assert.AreEqual(num0, at.Alpha0_dTau2, 1e-6 * Math.Max(1.0, Math.Abs(num0)), "Alpha0_dTau2");
        Assert.AreEqual(numR, at.AlphaR_dTau2, 1e-6 * Math.Max(1.0, Math.Abs(numR)), "AlphaR_dTau2");
    }

    [DataTestMethod]
    [DynamicData(nameof(TauDerivativePoints), DynamicDataSourceType.Method)]
    public void MixedDerivative_MatchesFiniteDifference_InBothDirections(double T, double delta)
    {
        double tau = Tc / T;
        double hTau = 1e-6 * tau;
        double hDelta = AdaptiveEps(delta);

        // d/dtau of alphaR_delta and d/ddelta of alphaR_tau must both equal alphaR_deltatau
        double viaTau = (StateAtTau(tau + hTau, delta).AlphaR_dDelta - StateAtTau(tau - hTau, delta).AlphaR_dDelta) / (2.0 * hTau);
        double viaDelta = (StateAtTau(tau, delta + hDelta).AlphaR_dTau - StateAtTau(tau, delta - hDelta).AlphaR_dTau) / (2.0 * hDelta);
        double analytic = StateAtTau(tau, delta).AlphaR_dDeltadTau;

        Assert.AreEqual(viaTau, analytic, 1e-6 * Math.Max(1.0, Math.Abs(analytic)), "via tau");
        Assert.AreEqual(viaDelta, analytic, 1e-6 * Math.Max(1.0, Math.Abs(analytic)), "via delta");
    }
}
