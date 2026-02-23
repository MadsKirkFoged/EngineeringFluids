using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace EngineeringFluids.Helmholtz
{
    public static class Saturation
    {
        public static readonly double Tr = 405.56;
        public static readonly double Tmax = 405.56;
        public static readonly double Tmin = 195.495;
        public static readonly double reducingValue = 11365000; // pc in Pa
        public static readonly double maxAbsErrorPercentage = 0.05181632089212851;

        public static readonly List<float> Coefficients =
        [
            -7.2257f,
            1.4263f,
            -0.59642f,
            -2.798f,
            -3.7869f
        ];

        public static readonly List<float> Exponents =
        [
            1.0f,
            1.5f,
            2.0f,
            3.6f,
            15.5f
        ];

        public static Pressure CalculateSaturationPressure(Temperature temperature)
        {

           double temperatureDouble = temperature.Kelvin;

            double theta = 1.0 - temperatureDouble / Tr;
            double sum = 0.0;

            for (int i = 0; i < Coefficients.Count; i++)
            {
                sum += Coefficients[i] * Math.Pow(theta, Exponents[i]);
            }


            double result = reducingValue * Math.Exp(Tr / temperatureDouble * sum);
            return Pressure.FromSI(result);

        }

        public static double CalculateSaturationPressureDouble(double temperature)
        {

            double temperatureDouble = temperature;

            double theta = 1.0 - temperatureDouble / Tr;
            double sum = 0.0;

            for (int i = 0; i < Coefficients.Count; i++)
            {
                sum += Coefficients[i] * Math.Pow(theta, Exponents[i]);
            }


            double result = reducingValue * Math.Exp(Tr / temperatureDouble * sum);
            return result;

        }

        public static float CalculateSaturationPressureFloat(float temperature)
        {

            float temperatureDouble = temperature;

            float theta = 1.0f - temperatureDouble / (float)Tr;
            float sum = 0.0f;

            for (int i = 0; i < Coefficients.Count; i++)
            {
                sum += Coefficients[i] * MathF.Pow(theta, Exponents[i]);
            }


            float result = (float)reducingValue * MathF.Exp((float)Tr / temperatureDouble * sum);
            return result;

        }





        public static Temperature CalculateSaturationTemperature(Pressure Pressure)
        {
            double low = Tmin;
            double high = Tmax;
            double epsilon = 0.0001; // Define your own tolerance level

            while (high - low > epsilon)
            {
                double mid = (low + high) / 2;
                Temperature midTemperature = Temperature.FromKelvin(mid);
                Pressure midPressure = CalculateSaturationPressure(midTemperature);

                if (Math.Abs(midPressure.SI - Pressure.SI) < epsilon)
                {
                    return midTemperature;
                }
                else if (midPressure.SI < Pressure.SI)
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
