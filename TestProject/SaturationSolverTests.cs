using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;

namespace TestProject;

[TestClass]
public class SaturationSolverTests
{
    [TestMethod]
    public void SaturationSolveAtT_EnforcesEquilibrium()
    {
        var T = Temperature.FromKelvin(400.0);

        var sat = new Ammonia().SolveAtT(T);

        // Build states at the solved densities
        var L = new Ammonia { Temperature = T, Density = sat.RhomolarL * new Ammonia().MolarMass };
        var V = new Ammonia { Temperature = T, Density = sat.RhomolarV * new Ammonia().MolarMass };

        // Check equilibrium: pressures match and lnphi matches
        Assert.AreEqual(L.Pressure.Pascal, V.Pressure.Pascal, 1e-2); // Pa tolerance
        Assert.AreEqual(L.LNFugacityCoefficient, V.LNFugacityCoefficient, 1e-10);

        // Also check ordering
        Assert.IsTrue(sat.RhomolarL.MolesPerCubicMeter > sat.RhomolarV.MolesPerCubicMeter);
    }
}