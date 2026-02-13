// See https://aka.ms/new-console-template for more information


// Given data
//using static System.Runtime.InteropServices.JavaScript.JSType;

using EngineeringFluids;
using EngineeringFluids.Fluids;
using EngineeringUnits;
using System.Diagnostics;
using EngineeringFluids.Helmholtz.Solvers;

var input = new Ammonia();
var inputD = new AmmoniaDouble();

for (int i = 0; i < 10000; i++)
{
    inputD.UpdatePT(Pressure.FromBar(21), Temperature.FromDegreeCelsius(100));
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



Debug.Print("");
