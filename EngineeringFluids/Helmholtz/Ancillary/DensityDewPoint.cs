using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineeringFluids.Helmholtz.Ancillary
{
    public static class VaporDensity
    {
        public static readonly double Tr = 405.56;
        public static readonly double Tmax = 405.56;
        public static readonly double Tmin = 195.495;
        public static readonly double reducingValue = 13696;
        public static readonly double maxAbsErrorPercentage = 0.312705;

        public static readonly List<double> Coefficients = new List<double>
        {
            -0.053296,
            -3.4589,
            -6.7572,
            -17.260,
            -43.120,
            -115.18
        };

        public static readonly List<double> Exponents = new List<double>
        {
            0.14,
            0.44,
            1.314,
            3.225,
            6.4,
            14.0
        };

        public static Molarity CalculateDensity(Temperature temperature)
        {

            if (temperature is null)
                return null;

            var temperatureAsDouble = temperature.Kelvin;

            double THETA = 1.0 - temperatureAsDouble / Tr;
            double summer = 0.0;

            for (int i = 0; i < Coefficients.Count; i++)
            {
                summer += Coefficients[i] * Math.Pow(THETA, Exponents[i]);
            }

            //return reducingValue * Math.Exp(summer);

            double result = reducingValue * Math.Exp(summer);
            return Molarity.FromMolesPerCubicMeter(result);

        }

        public static double CalculateDensityDouble(double temperature)
        {

            var temperatureAsDouble = temperature;

            double THETA = 1.0 - temperatureAsDouble / Tr;
            double summer = 0.0;

            for (int i = 0; i < Coefficients.Count; i++)
            {
                summer += Coefficients[i] * Math.Pow(THETA, Exponents[i]);
            }

            //return reducingValue * Math.Exp(summer);

            double result = reducingValue * Math.Exp(summer);
            return result;

        }

        public static Temperature CalculateTemperature(Molarity MolarDensity)
        {

            if (MolarDensity is null)
                return null;

            double low = Tmin;
            double high = Tmax;
            double epsilon = 0.0001; // Define your own tolerance level

            while (high - low > epsilon)
            {
                double mid = (low + high) / 2;
                Temperature midTemperature = Temperature.FromKelvin(mid);
                Molarity midDensity = CalculateDensity(midTemperature);

                if (Math.Abs(midDensity.MolesPerCubicMeter - MolarDensity.MolesPerCubicMeter) < epsilon)
                {
                    return midTemperature;
                }
                else if (midDensity.MolesPerCubicMeter < MolarDensity.MolesPerCubicMeter)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }

            return Temperature.FromKelvin((low + high) / 2);
        }
    }


}
