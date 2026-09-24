using EngineeringFluids.Fluids;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace TestProject;

[TestClass]
public class DerivativeTests
{
    private static readonly Ammonia Ref = new();
    private static readonly double RhoCritMass = Ref.Critical.MolarDensity * Ref.MolarMass; // kg/m3

    private static double AdaptiveEps(double delta)
    {
        double eps = 1e-6 * Math.Max(1.0, Math.Abs(delta));
        if (delta - eps <= 0.0)
            eps = 0.5 * delta; // keep delta-eps > 0
        return eps;
    }

    private static Ammonia StateAt(double T, double delta)
        => new Ammonia { Temperature = T, Density = RhoCritMass * delta };

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
        AssertDDelta2MatchesFiniteDifference(400.0, 9.0 / RhoCritMass, "Vapor-like point");
    }

    [TestMethod]
    public void AlphaR_dDelta2_MatchesFiniteDifference_LiquidLikePoint()
    {
        AssertDDelta2MatchesFiniteDifference(400.0, 1.4748450193131541, "Liquid-like point");
    }
}
