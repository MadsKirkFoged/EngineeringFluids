using System;

namespace EngineeringFluids.Helmholtz;
public static class IdealGasHelmholtzLead
{
    // alpha0 = ln(delta) + a1 + a2*tau
    public static IdealDerivatives Derivatives(double delta, double tau)
    {
        double value = Math.Log(delta) + -6.59406093943886 + 5.60101151987913 * tau;
        return new IdealDerivatives(value, 5.60101151987913, 0.0);
    }
}
