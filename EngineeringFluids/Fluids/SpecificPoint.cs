using EngineeringUnits;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineeringFluids.Fluids;
public  record SpecificPoint
{
    public required Temperature Temperature { get; init; }

    public required double MolarEnthalpy { get; init; }

    public required Pressure Pressure { get; init; }

    public required Molarity MolarDensity { get; init; }

    public required MolarEntropy MolarEntropy { get; init; }

}
