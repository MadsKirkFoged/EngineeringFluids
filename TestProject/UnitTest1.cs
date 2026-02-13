using EngineeringFluids.Fluids;
using EngineeringFluids.Helmholtz;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringUnits;
using EngineeringUnits.Units;
using static EngineeringFluids.Helmholtz.Phase;

namespace TestProject;

[TestClass]
public class AmmoniaTests
{
    //[TestMethod]
    //public void ComparingProperties()
    //{

    //    //Arrange
    //    var test = new Ammonia();


    //    //Act
    //    test.Temperature =  Temperature.FromKelvin(400.0);
    //    test.Density = Density.FromKilogramPerCubicMeter(9);


    //    //Assert
    //    Assert.AreEqual(1.0139, test.Tau);
    //    Assert.AreEqual(0.03858520868560085, test.Delta, 0.00000000000000001);
    //    Assert.AreEqual(-4.1654333380781647, test.Alpha0);
    //    Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.00000000000000010);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau);
    //    Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.00000000000001000);
    //    Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.00000000000000010);

    //    Assert.AreEqual(1654042.55202634, test.Pressure.SI, 0.00000001);
    //    Assert.AreEqual(6298.33191543816, test.Entropy.SI, 0.0000000001);
    //    Assert.AreEqual(1694261.00050462, test.InternalEnergy.SI, 0.00000001);
    //    Assert.AreEqual(1878043.50628532, test.Enthalpy.SI, 0.00000001);

    //    Assert.AreEqual(31984.057494662303592976800m, test.HMolarEnthalpy.AsSI, 0.00000001m);
    //    Assert.AreEqual(23095.272807360314, test.MolarEnthalpyDewPointGuess, 0.00000001);
    //    Assert.AreEqual(18029.539580378132, test.MolarEnthalpyBubblePointGuess, 0.00000001);

    //    Assert.AreEqual(1088593.6022863, test.EnthalpyBubblePoint.SI, 10000);
    //    Assert.AreEqual(1434784.01069596, test.EnthalpyDewPoint.SI, 10000);

    //    Assert.AreEqual(Phases.Gas, test.Phase);

    //    Assert.AreEqual(-1, test.Quality);

    //    Assert.AreEqual(528.46301815798927, test.MolarDensity.MolesPerCubicMeter, 0.00000001);


    //    Assert.AreEqual(20199.47738451296, test.MolarDensityBubblePoint.SI, 100);
    //    Assert.AreEqual(7685.3662399232844, test.MolarDensityDewPoint.SI, 100);

    //    Assert.AreEqual(10297199.3950676, test.SatPressure.SI, 301693);

    //    //-10921.489566340852
    //    //-10921.48956634085

    //}




    [TestMethod]
    public void UpdatePT()
    {

        //Arrange
        var test = new Ammonia();
        Pressure p = Pressure.FromBar(16.5404255202634);
        Temperature t = Temperature.FromKelvin(400);


        //Act
        test.UpdatePT(p, t);

        var d2 = test.AlphaR_dDelta2;

        //Assert
        Assert.AreEqual(1.0139, test.Tau);
        Assert.AreEqual(0.03858520868560085, test.Delta, 0.000001);
        Assert.AreEqual(-4.1654333380781647, test.Alpha0, 0.000001);
        Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.000001);

        Assert.AreEqual(8.71121706208864, test.Alpha0_dTau);
        Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.0000001000);
        Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.0000001);

        Assert.AreEqual(1654042.55202634, test.Pressure.SI, 1);
        Assert.AreEqual(6298.33191543816, test.Entropy.SI, 0.001);
        Assert.AreEqual(1694261.00050462, test.InternalEnergy.SI, 0.01);
        Assert.AreEqual(1878043.50628532, test.Enthalpy.SI, 0.1);


        Assert.AreEqual(Phases.Gas, test.Phase);
        Assert.AreEqual(-1, test.Quality);



    }

    //[TestMethod]
    //public void UpdateDP()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    Pressure p = Pressure.FromBar(16.5404255202634);
    //    Density d = test.Density = Density.FromKilogramPerCubicMeter(9);


    //    //Act
    //    test.UpdateDP(d,p);


    //    //Assert
    //    Assert.AreEqual(400, test.Temperature.SI, 0.00000001);
    //    Assert.AreEqual(1.0139, test.Tau, 0.00000001);
    //    Assert.AreEqual(0.03858520868560085, test.Delta, 0.00000000000001);
    //    Assert.AreEqual(-4.1654333380781647, test.Alpha0, 0.0000000001);
    //    Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.000000000001);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau, 0.0000000001);
    //    Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.000000001);
    //    Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.0000000001);

    //    Assert.AreEqual(1654042.55202634, test.Pressure.SI, 0.0001);
    //    Assert.AreEqual(6298.33191543816, test.Entropy.SI, 0.000001);
    //    Assert.AreEqual(1694261.00050462, test.InternalEnergy.SI, 0.00001);
    //    Assert.AreEqual(1878043.50628532, test.Enthalpy.SI, 0.0001);


    //    Assert.AreEqual(Phases.Gas, test.Phase);
    //    Assert.AreEqual(-1, test.Quality);



    //}

    //[TestMethod]
    //public void UpdateTQ()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    Temperature t = Temperature.FromKelvin(400);
    //    double q = 1;


    //    //Act
    //    test.UpdateTQ(t, q);


    //    //Assert
    //    Assert.AreEqual(400, test.Temperature.SI, 0.00000001);
    //    Assert.AreEqual(1.0139, test.Tau, 0.00000001);
    //    Assert.AreEqual(0.03858520868560085, test.Delta, 0.00000000000001);
    //    Assert.AreEqual(-4.1654333380781647, test.Alpha0, 0.0000000001);
    //    Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.000000000001);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau, 0.0000000001);
    //    Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.000000001);
    //    Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.0000000001);

    //    Assert.AreEqual(1654042.55202634, test.Pressure.SI, 0.0001);
    //    Assert.AreEqual(6298.33191543816, test.Entropy.SI, 0.000001);
    //    Assert.AreEqual(1694261.00050462, test.InternalEnergy.SI, 0.00001);
    //    Assert.AreEqual(1878043.50628532, test.Enthalpy.SI, 0.0001);


    //    Assert.AreEqual(Phases.Gas, test.Phase);
    //    Assert.AreEqual(-1, test.Quality);



    //}

    //[TestMethod]
    //public void UpdatePQ()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    Pressure p = Pressure.FromPascal(10297199.3950676);
    //    double q = 0;


    //    //Act
    //    test.UpdatePQ(p, q);


    //    //Assert
    //    Assert.AreEqual(400, test.Temperature.SI, 0.01);
    //    Assert.AreEqual(1.0139, test.Tau, 0.001);
    //    Assert.AreEqual(1.4748450193131541, test.Delta, 0.01);
    //    Assert.AreEqual(-0.52199415483614919, test.Alpha0, 0.01);
    //    Assert.AreEqual(-1.4502590354825864, test.AlphaR, 0.01);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau, 0.1);
    //    Assert.AreEqual(-0.5741079073022356, test.AlphaR_dDelta, 0.01);
    //    Assert.AreEqual(-3.3643980447855943, test.AlphaR_dTau, 0.01);

    //    Assert.AreEqual(10297199.3950676, test.Pressure.SI, 0.0001);
    //    Assert.AreEqual(3609.5238910657, test.Entropy.SI, 1);
    //    Assert.AreEqual(1058660.5447384, test.InternalEnergy.SI, 300);
    //    Assert.AreEqual(1088593.6022863, test.Enthalpy.SI, 300);


    //    //Assert.AreEqual(Phases.Liquid, test.Phase);
    //    Assert.AreEqual(0, test.Quality);



    //}

    //[TestMethod]
    //public void AncillaryEnthalpyPoints()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    Temperature t = Temperature.FromKelvin(405);


    //    //Act
    //    test.Temperature = t;


    //    //Assert
    //    Assert.AreEqual(8009.2448070352784, test.MolarEnthalpyBubblePointGuess, 10000);
    //    Assert.AreEqual(25569.570962501886, test.MolarEnthalpyDewPointGuess, 10000);



    //}

    //[TestMethod]
    //public void Fugacity()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    Temperature t = Temperature.FromDegreeCelsius(26);
    //    //Pressure p = Pressure.FromBar(10.6);
    //    Density d = Density.FromKilogramPerCubicMeter(601.5);

    //    //Act
    //    test.Temperature = t; 
    //    test.Density = d;



    //    //Assert

    //    Assert.AreEqual(0.012071562946959769, test.CompressibilityFactor, 0);
    //    Assert.AreEqual(-0.7149694222776300, test.Alpha, 0);
    //    Assert.AreEqual(-3.5421648633914908, test.AlphaR, 0);
    //    Assert.AreEqual(2.5786929629766386, test.Delta, 0);
    //    Assert.AreEqual(-0.38311396066125958, test.AlphaR_dDelta, 0);
    //    Assert.AreEqual(1.3557078388768178, test.Tau, 0);
    //    Assert.AreEqual(-5.4962222891656625, test.AlphaR_dTau, 0);


    //    Assert.AreEqual(0.87281238244676373, test.FugacityCoefficient, 0);
    //    Assert.AreEqual(0.87281238244676373, test.LNFugacityCoefficient, 0);

    //    Assert.AreEqual(925181.12539376516, test.Fugacity.Pascal, 0);


    //}

    [TestMethod]
    public void LiqAndGasFinder()
    {

        //Arrange
        var test = new Ammonia();
        Temperature t = Temperature.FromDegreeCelsius(26);
        Density d = Density.FromKilogramPerCubicMeter(550);

        //Act
        test.Temperature = t;
        test.Density = d;

        test.FindEnthalpyAtLiqAndGas();

        //Assert

        //Assert.AreEqual(0.012071562946959769, test.CompressibilityFactor, 0);
        //Assert.AreEqual(-0.7149694222776300, test.Alpha, 0);
        //Assert.AreEqual(-3.5421648633914908, test.AlphaR, 0);
        //Assert.AreEqual(2.5786929629766386, test.Delta, 0);
        //Assert.AreEqual(-0.38311396066125958, test.AlphaR_dDelta, 0);
        //Assert.AreEqual(1.3557078388768178, test.Tau, 0);
        //Assert.AreEqual(-5.4962222891656625, test.AlphaR_dTau, 0);


        //Assert.AreEqual(0.87281238244676373, test.FugacityCoefficient, 0);
        //Assert.AreEqual(0.87281238244676373, test.LNFugacityCoefficient, 0);

        //Assert.AreEqual(925181.12539376516, test.Fugacity.Pascal, 0);


    }


    [TestMethod]
    public void UpdatePT_Newton_VaporPoint()
    {
        var test = new Ammonia();
        Pressure p = Pressure.FromBar(16.5404255202634);
        Temperature t = Temperature.FromKelvin(400);

        test.UpdatePT(p, t);

        Assert.AreEqual(400, test.Temperature.Kelvin, 1e-10);
        Assert.AreEqual(1654042.55202634, test.Pressure.Pascal, 1);
        Assert.AreEqual(0.03858520868560085, test.Delta, 0.00001);
    }

    //[TestMethod]
    //public void UpdatePQ_UsesEOSDome()
    //{
    //    var a = new Ammonia();
    //    var p = Pressure.FromPascal(10297199.3950676); // your earlier value
    //    a.UpdatePQ(p, 0.0);

    //    Assert.AreEqual(0.0, a.Quality, 0.0);
    //    Assert.AreEqual(p.Pascal, a.Pressure.Pascal, 1e-1);

    //    // For q=0, properties should match saturated liquid endpoint exactly
    //    var sat = new Ammonia().SolveAtP(p);
    //    var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * new Ammonia().MolarMass };

    //    Assert.AreEqual(L.Enthalpy.SI, a.Enthalpy.SI, 1e-6);
    //}

    [TestMethod]
    public void SolveAtP_DoesNotTouchTriplePoint_ForHighPressure()
    {
        var P = Pressure.FromPascal(10297199.3950676); // ~10 MPa
        var sat = new Ammonia().SolveAtP(P);

        Assert.IsTrue(sat.T.Kelvin > 300.0);
        Assert.IsTrue(sat.T.Kelvin < new Ammonia().Critical.Temperature.Kelvin);
        Assert.AreEqual(P.Pascal, sat.Psat.Pascal, 5.0); // loose until everything is tuned
    }

    //[TestMethod]
    //public void UpdatePH_InsideDome_ReturnsSameEnthalpy()
    //{
    //    Pressure p = Pressure.FromBar(102.97); //I like bar better (What I normally work with)

    //    var sat = new Ammonia().SolveAtP(p);

    //    var L = new Ammonia { Temperature = sat.T, Density = sat.RhomolarL * new Ammonia().MolarMass };
    //    var V = new Ammonia { Temperature = sat.T, Density = sat.RhomolarV * new Ammonia().MolarMass };

    //    double q = 0.3;
    //    var hTarget = (1 - q) * L.Enthalpy + q * V.Enthalpy;

    //    var a = new Ammonia();
    //    a.UpdatePH(p, hTarget);

    //    Assert.AreEqual(sat.T.Kelvin, a.Temperature.Kelvin, 1e-3);
    //    Assert.AreEqual(p.Pascal, a.Pressure.Pascal, 1e-1);
    //    Assert.AreEqual(hTarget.AsSI, a.Enthalpy.AsSI, 1e-2m);
    //}

    //[TestMethod]
    //public void UpdatePH_SinglePhase_ReconstructsKnownState()
    //{
    //    // Known state from your earlier test
    //    var T = Temperature.FromKelvin(400);
    //    var rho = Density.FromKilogramPerCubicMeter(9);

    //    var refState = new Ammonia { Temperature = T, Density = rho };
    //    var p = refState.Pressure;
    //    var h = refState.Enthalpy;

    //    var a = new Ammonia();
    //    a.UpdatePH(p, h);

    //    Assert.AreEqual(T.Kelvin, a.Temperature.Kelvin, 1e-2);
    //    Assert.AreEqual(p.Pascal, a.Pressure.Pascal, 1);
    //    Assert.AreEqual(h.SI, a.Enthalpy.SI, 1e-1);
    //}




    //[TestMethod]
    //public void UpdateDQGas()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    double q = 1;
    //    Density d = test.Density = Density.FromKilogramPerCubicMeter(9);


    //    //Act
    //    test.UpdateDQ(d, q);


    //    //Assert
    //    Assert.AreEqual(400, test.Temperature.SI, 0.00000001);
    //    Assert.AreEqual(1.0139, test.Tau, 0.00000001);
    //    Assert.AreEqual(0.03858520868560085, test.Delta, 0.00000000000001);
    //    Assert.AreEqual(-4.1654333380781647, test.Alpha0, 0.0000000001);
    //    Assert.AreEqual(-0.059555880911521322, test.AlphaR, 0.000000000001);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau, 0.0000000001);
    //    Assert.AreEqual(-1.5263639324940277, test.AlphaR_dDelta, 0.000000001);
    //    Assert.AreEqual(-0.15426556068786079, test.AlphaR_dTau, 0.0000000001);

    //    Assert.AreEqual(1654042.55202634, test.Pressure.SI, 0.0001);
    //    Assert.AreEqual(6298.33191543816, test.Entropy.SI, 0.000001);
    //    Assert.AreEqual(1694261.00050462, test.InternalEnergy.SI, 0.00001);
    //    Assert.AreEqual(1878043.50628532, test.Enthalpy.SI, 0.0001);


    //    Assert.AreEqual(Phases.Gas, test.Phase);
    //    Assert.AreEqual(1, test.Quality);



    //}

    //[TestMethod]
    //public void UpdateDQLiq()
    //{

    //    //Arrange
    //    var test = new Ammonia();
    //    double q = 0;
    //    Density d = test.Density = Density.FromKilogramPerCubicMeter(9);


    //    //Act
    //    test.UpdateDQ(d, q);


    //    //Assert
    //    Assert.AreEqual(400, test.Temperature.SI, 0.00000001);
    //    Assert.AreEqual(1.0139, test.Tau, 0.00000001);
    //    Assert.AreEqual(1.4748450193131541, test.Delta, 0.00000000000001);
    //    Assert.AreEqual(-0.52199415483614919, test.Alpha0, 0.0000000001);
    //    Assert.AreEqual(-1.4502590354825864, test.AlphaR, 0.000000000001);

    //    Assert.AreEqual(8.71121706208864, test.Alpha0_dTau, 0.0000000001);
    //    Assert.AreEqual(-0.5741079073022356, test.AlphaR_dDelta, 0.000000001);
    //    Assert.AreEqual(-3.3643980447855943, test.AlphaR_dTau, 0.0000000001);

    //    Assert.AreEqual(10297199.3950676, test.Pressure.SI, 0.0001);
    //    Assert.AreEqual(3609.5238910657, test.Entropy.SI, 0.000001);
    //    Assert.AreEqual(1058660.5447384, test.InternalEnergy.SI, 0.00001);
    //    Assert.AreEqual(1088593.6022863, test.Enthalpy.SI, 0.0001);


    //    Assert.AreEqual(Phases.Liquid, test.Phase);
    //    Assert.AreEqual(0, test.Quality);



    //}
}
