using EngineeringFluids.Fluids;
using EngineeringUnits;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace TestProject;

[TestClass]
public class DerivativeTests
{
    private static double AdaptiveEps(double delta)
    {
        double eps = 1e-6 * Math.Max(1.0, Math.Abs(delta));
        if (delta - eps <= 0.0)
        {
            eps = 0.5 * delta; // keep delta-eps > 0
        }
        return eps;
    }

    private static double AlphaR_dDelta_at(Temperature T, double delta)
    {
        // delta = rho_mass / rho_crit_mass, so rho_mass = delta * rho_crit_mass
        var rhoCritMass = new Ammonia().Critical.MolarDensity * new Ammonia().MolarMass; // Density
        var a = new Ammonia
        {
            Temperature = T,
            Density = rhoCritMass * delta
        };
        return a.AlphaR_dDelta;
    }

    private static double AlphaR_dDelta2_at(Temperature T, double delta)
    {
        var rhoCritMass = new Ammonia().Critical.MolarDensity * new Ammonia().MolarMass; // Density
        var a = new Ammonia
        {
            Temperature = T,
            Density = rhoCritMass * delta
        };
        return a.AlphaR_dDelta2;
    }

    [TestMethod]
    public void AlphaR_dDelta2_MatchesFiniteDifference_VaporLikePoint()
    {
        // Use an existing vapor-like point from your own tests:
        // T = 400 K, rho_mass = 9 kg/m3 => delta ~ 0.0386 in your assertions
        Temperature T = Temperature.FromKelvin(400.0);

        // Compute delta from the same definition your Ammonia uses
        var rhoCritMass = new Ammonia().Critical.MolarDensity * new Ammonia().MolarMass; // Density
        double delta = (double)(Density.FromKilogramPerCubicMeter(9.0) / rhoCritMass).AsSI;

        double eps = AdaptiveEps(delta);

        double num = (AlphaR_dDelta_at(T, delta + eps) - AlphaR_dDelta_at(T, delta - eps)) / (2.0 * eps);
        double ana = AlphaR_dDelta2_at(T, delta);

        // Tolerance: start moderately loose; tighten once you see it’s stable
        // If this fails badly, it usually means a sign error in one term.
        Assert.AreEqual(num, ana, 1e-7, $"Vapor-like point mismatch: num={num}, ana={ana}, delta={delta}, eps={eps}");
    }

    [TestMethod]
    public void AlphaR_dDelta2_MatchesFiniteDifference_LiquidLikePoint()
    {
        // Use a liquid-like delta around 1-3; your UpdatePQ test shows delta ~ 1.47 at ~400 K
        Temperature T = Temperature.FromKelvin(400.0);
        double delta = 1.4748450193131541;

        double eps = AdaptiveEps(delta);

        double num = (AlphaR_dDelta_at(T, delta + eps) - AlphaR_dDelta_at(T, delta - eps)) / (2.0 * eps);
        double ana = AlphaR_dDelta2_at(T, delta);

        Assert.AreEqual(num, ana, 1e-7, $"Liquid-like point mismatch: num={num}, ana={ana}, delta={delta}, eps={eps}");
    }
}