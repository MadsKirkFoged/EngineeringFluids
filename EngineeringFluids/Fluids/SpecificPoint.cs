using EngineeringUnits.Fast;

namespace EngineeringFluids.Fluids;

public record SpecificPoint
{
    public required Temperature Temperature { get; init; }

    public required MolarEnergy MolarEnthalpy { get; init; }

    public required Pressure Pressure { get; init; }

    public required Molarity MolarDensity { get; init; }

    public required MolarEntropy MolarEntropy { get; init; }
}
