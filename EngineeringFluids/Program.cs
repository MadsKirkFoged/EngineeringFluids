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





