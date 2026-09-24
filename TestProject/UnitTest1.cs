using EngineeringFluids.Fluids;
using static EngineeringFluids.Helmholtz.Phase;

namespace TestProject;

[TestClass]
public class AmmoniaTests
{
    private const double P = 1654042.55202634; // Pa (16.5404255202634 bar)
    private const double T = 400.0;             // K

    [TestMethod]
    public void UpdatePT()
    {
        var test = new Ammonia();

        test.UpdatePT(P, T);

        Assert.AreEqual(1.0139, test.Tau);
        Assert.AreEqual(0.03858520868560085, test.Delta, 0.000001);
        Assert.AreEqual(-4.1654333380781647, test.Alpha0, 0.000001);
        Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.000001);

        Assert.AreEqual(8.71121706208864, test.Alpha0_dTau);
        Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.0000001000);
        Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.0000001);

        Assert.AreEqual(1654042.55202634, test.Pressure, 1);
        Assert.AreEqual(6298.33191543816, test.Entropy, 0.001);
        Assert.AreEqual(1694261.00050462, test.InternalEnergy, 0.01);
        Assert.AreEqual(1878043.50628532, test.Enthalpy, 0.1);

        Assert.AreEqual(Phases.Gas, test.Phase);
        Assert.AreEqual(-1, test.Quality);
    }

    [TestMethod]
    public void UpdatePT_Newton_VaporPoint()
    {
        var test = new Ammonia();

        test.UpdatePT(P, T);

        Assert.AreEqual(400, test.Temperature, 1e-10);
        Assert.AreEqual(1654042.55202634, test.Pressure, 1);
        Assert.AreEqual(0.03858520868560085, test.Delta, 0.00001);
    }
}
