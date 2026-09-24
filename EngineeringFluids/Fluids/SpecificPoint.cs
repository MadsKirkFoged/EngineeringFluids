namespace EngineeringFluids.Fluids;

public record SpecificPoint
{
    public required double Temperature { get; init; }

    public required double MolarEnthalpy { get; init; }

    public required double Pressure { get; init; }

    public required double MolarDensity { get; init; }

    public required double MolarEntropy { get; init; }
}
