namespace EngineeringFluids.Helmholtz;

/// <summary>
/// An ideal-gas Helmholtz term and its tau derivatives at one (delta, tau). The ideal part has no delta
/// dependence beyond ln(delta), so these are the only derivatives the properties need.
/// </summary>
public readonly record struct IdealDerivatives(double Value, double dTau, double dTau2);

/// <summary>
/// A residual Helmholtz term and every derivative the properties need, evaluated together at one (delta, tau)
/// so the expensive exp/log factors are computed once per state instead of once per derivative.
/// </summary>
public readonly record struct ResidualDerivatives(
    double Value,
    double dDelta,
    double dTau,
    double dDelta2,
    double dTau2,
    double dDeltadTau);
