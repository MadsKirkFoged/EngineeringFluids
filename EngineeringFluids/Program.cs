// See https://aka.ms/new-console-template for more information


// Given data
//using static System.Runtime.InteropServices.JavaScript.JSType;

using EngineeringFluids;
using EngineeringFluids.Fluids;
using EngineeringUnits;
using System.Diagnostics;
using EngineeringFluids.Helmholtz.Solvers;
using EngineeringFluids.Helmholtz;
using SharpFluids;
using EngineeringFluids.Helmholtz.Ancillary;
using System.Globalization;






//var input = new Ammonia();
var inputD = new AmmoniaDouble();


double pressure = 2100000; // Pressure.FromBar(21).SI;
double temperature = 373.15; // Temperature.FromDegreeCelsius(100).SI;

inputD.UpdatePT(pressure, temperature);

double T = 400.0; // Temperature in K
 double rho = 9.0; // Density in kg/m^3

 double Tc = 405.56; // Critical temperature in K
 double rhoc = 13696.0; // Critical density in mol/m^3
 double M = 0.01703056; // Molar mass in kg/mol
 double R = 8.3144598; // Gas constant in J/(mol*K)

 double tau = Tc / T;
 double delta = rho / (rhoc * M);


var check1 = ResidualHelmholtzGaussian.alphaR_dDelta2(delta, tau);
var check2 = ResidualHelmholtzGaussianFast.alphaR_dDelta2(delta, tau);
var check21 = ResidualHelmholtzGaussianSIMD.alphaR_dDelta2(delta, tau);

var check3 = ResidualHelmholtzGaussian.alphaR_dDelta(delta, tau);
var check4 = ResidualHelmholtzGaussianFast.alphaR_dDelta(delta, tau);
var check41 = ResidualHelmholtzGaussianSIMD.alphaR_dDelta(delta, tau);

var check5 = ResidualHelmholtzGaussian.alphaR_dTau(delta, tau);
var check6 = ResidualHelmholtzGaussianFast.alphaR_dTau(delta, tau);
var check61 = ResidualHelmholtzGaussianSIMD.alphaR_dTau(delta, tau);

var check7 = ResidualHelmholtzGaussian.alphaR2_dTau(delta, tau);
var check8 = ResidualHelmholtzGaussianFast.alphaR2_dTau(delta, tau);
var check81 = ResidualHelmholtzGaussianSIMD.alphaR2_dTau(delta, tau);

var check9 = ResidualHelmholtzGaussian.alphaR_dDelta2(delta, tau);
var check10 = ResidualHelmholtzGaussianFast.alphaR_dDelta2(delta, tau);
var check101 = ResidualHelmholtzGaussianSIMD.alphaR_dDelta2(delta, tau);

//Fluid Bubble = new Fluid(FluidList.Ammonia);
//Fluid Dew = new Fluid(FluidList.Ammonia);

//var test = FastPressurePoly.Pressure(373.15f);
//var test2 = Saturation.CalculateSaturationPressureDouble(373.15d);
//Bubble.UpdateXT(0, Temperature.FromKelvin(373.15));




//for (double i = 195.495; i < 405.56; i = i + 0.005)
//{
    //Bubble.UpdateXT(0, Temperature.FromKelvin(i));
    //Dew.UpdateXT(1, Temperature.FromKelvin(i));

    
    //var DewEnthalpyAncillary = AncillaryEnthalpyDewPoint.MolarEnthalpy(Temperature.FromSI(i));

    //var BubbleDensityAncillary = LiquidDensity.CalculateDensityDouble(i);
    ////var DewDensityAncillary = VaporDensity.CalculateDensityDouble(i);


    ////var BubbleEnthalpy = Bubble.Enthalpy.SI * inputD.MolarMass;
    ////var BubbleEnthalpyfast = BubbleEnthalpyFast.BubbleEnthalpy((float)i);

    //var BD = Bubble.Density.SI / inputD.MolarMass;
    //var BDF = BubbleDensityFast.Density((float)i);

    //var newdiss = Math.Abs(BD - BDF);
    //Debug.Print(FormattableString.Invariant($"{i}  NewDiff: {newdiss}"));


    //var BubbleEnthalpyAncillary = AncillaryEnthalpyBubblePoint.MolarEnthalpy(Temperature.FromSI(i));

     

    //var BubbleDensity = Bubble.Density.SI / inputD.MolarMass;

    //var DewDensity = Dew.Density.SI / inputD.MolarMass;

    //var temperature = i;
    //var pressure = Bubble.Pressure.SI;



    //Write(
    //    filePath: "output/data2.csv",
    //    temperature: temperature,
    //    pressure: pressure,
    //    BubbleEnthalpy: BubbleEnthalpy,
    //    DewEnthalpy: DewEnthalpy,
    //    BubbleDensity: BubbleDensity,
    //    DewDensity: DewDensity,
    //    append: true  // set to false to overwrite
    //);




    //var diffPressure = myfluid2.Pressure - myfluid.Pressure;
    //var diffTemp = myfluid3.Temperature - myfluid4.Temperature;

    //var diffRoundtripTemp = myfluid3.Temperature - myfluid.Temperature;
    //var diffRoundtripPres = myfluid3.Pressure - myfluid.Pressure;

    //Debug.Print(FormattableString.Invariant($"1: {diffPressure?.AsSI}  2: {diffTemp?.AsSI}  3: {diffRoundtripTemp?.AsSI}  4: {diffRoundtripPres?.AsSI}"));

    //var olddiff = Math.Abs(Saturation.CalculateSaturationPressureDouble(i) - Bubble.Pressure.SI);
    //var newdiss = Math.Abs(FastPressurePoly.Pressure((float)i) - Bubble.Pressure.SI);
    //var newdiss2 = Math.Abs(SaturationPressureFast.Pressure((float)i) - Bubble.Pressure.SI);



    //Debug.Print(FormattableString.Invariant($"{i}  NewDiff: {newdiss}  NewDiff2: {newdiss2}"));

    //Debug.Print(FormattableString.Invariant($"{i};{myfluid.Pressure.AsSI}"));
//}



/// <summary>
/// Writes the values to a single line in a file using semicolons as separators
/// and CultureInfo.InvariantCulture for numeric formatting.
/// </summary>
static void Write(
    string filePath,
    double temperature,
    double pressure,
    double BubbleEnthalpy,
    double DewEnthalpy,
    double BubbleDensity,
    double DewDensity,
    bool append = false)
{
    var values = new[]
    {
            temperature,
            pressure,
            BubbleEnthalpy,
            DewEnthalpy,
            BubbleDensity,
            DewDensity
        };

    var line = string.Join(";", values.Select(v => v.ToString(CultureInfo.InvariantCulture)));

    // Create directory if needed
    var dir = Path.GetDirectoryName(filePath);
    if (!string.IsNullOrEmpty(dir) && !Directory.Exists(dir))
        Directory.CreateDirectory(dir);
    if (append)
        File.AppendAllText(filePath, line + Environment.NewLine);
    else
        File.WriteAllText(filePath, line + Environment.NewLine);
}






//Pressure p = Pressure.Parse("10.297 MPa");
//var result = input.SolveAtP(p);


//var span = (input.Critical.Temperature - input.TripleLiquid.Temperature);
//var step = span / 1000;

//var test = new Ammonia();





//for (int i = 0; i < 1000; i++)
//{
//    test.Temperature = input.Critical.Temperature - step*i;
//    Console.WriteLine($"{test.MolarEnthalpyDewPointGuess}");
//}






//double T = 400.0; // Temperature in K
//double rho = 9.0; // Density in kg/m^3


//double Tc = 405.56; // Critical temperature in K
//double rhoc = 13696.0; // Critical density in mol/m^3

//double M = 0.01703052; // Molar mass in kg/mol
//double R = 8.3144598; // Gas constant in J/(mol*K)

//// Calculate reduced temperature and density
//double tau = (Tc / T);                                          //Should be: 1.0139



//double delta = (rho / (rhoc * M));                              //Should be: 0.03858520868560085


//var alpha0NEW = Ammonia.alpha0(delta, tau);                  //Should be: -4.1654333380781647



//var alphaRNEW = Ammonia.alphaR(delta, tau);                  //should be: -0.059555880911521322



//var alpha0_dTauNEW = Ammonia.alpha0_dTau(delta, tau);        //Should be: 8.71121706208864



//var alphaR_dDeltaNEW = Ammonia.alphaR_dDelta(delta, tau);    //Should be: -1.5263639324940277


//var alphaR_dTauNEW = Ammonia.alphaR_dTau(delta, tau);        //should be: -0.15426556068786079

//ar test = 10;


//var alpha_dTau = alpha0_dTauNEW + alphaR_dTauNEW;

//double rhomolar = (rho / M);
//double pressure = (rhomolar * R * T * (1 + delta * alphaR_dDeltaNEW));




//////entropy

//var molarentropy = R * (tau * (alpha0_dTauNEW + alphaR_dTauNEW) - alpha0NEW - alphaRNEW);
//var entropy = molarentropy / M;



//////molar Internal Energy
//////R_u* T * tau * (da0_dTau + dar_dTau);

//var molarInternalEnergy = R * T * tau * (alpha0_dTauNEW + alphaR_dTauNEW);

//var InternalEnergy = molarInternalEnergy / M;


////// 1694261.00050462 
////// 1682328.8078847416


////// R_u * T * (1 + tau * (da0_dTau + dar_dTau) + delta * dar_dDelta);

//var hmolarmolarEnthalpy = R * T * (1 + tau * (alpha0_dTauNEW + alphaR_dTauNEW) + delta * alphaR_dDeltaNEW);
//var Enthalpy = hmolarmolarEnthalpy / M;


//1878043.50628532
//1866110.9434088822


//1654038.88
//1654042.55202634

//Skal vi kigge på at lave et Pressure-Temperature update?
// --> skal have kigget på hvordan coolprop gør det



//Vi skal have beregnet dalphar_dDelta 



//Debug.Print("");
