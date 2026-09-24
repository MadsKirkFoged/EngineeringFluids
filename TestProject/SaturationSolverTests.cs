using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz.Solvers;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace TestProject;

[TestClass]
public class SaturationSolverTests
{
    [TestMethod]
    public void SaturationSolveAtT_EnforcesEquilibrium()
    {
        const double T = 400.0;
        var reference = new Ammonia();

        var sat = reference.SolveAtT(T);

        // Build states at the solved densities
        var L = new Ammonia { Temperature = T, Density = sat.RhomolarL * reference.MolarMass };
        var V = new Ammonia { Temperature = T, Density = sat.RhomolarV * reference.MolarMass };

        // Check equilibrium: pressures match and lnphi matches
        Assert.AreEqual(L.Pressure, V.Pressure, 1e-2); // Pa tolerance
        Assert.AreEqual(L.LNFugacityCoefficient, V.LNFugacityCoefficient, 1e-10);

        Assert.IsTrue(sat.RhomolarL > sat.RhomolarV);
    }
}
