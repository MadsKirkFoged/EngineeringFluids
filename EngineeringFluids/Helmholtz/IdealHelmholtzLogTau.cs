using System;

namespace EngineeringFluids.Helmholtz;
public static class IdealHelmholtzLogTau
{
    // alpha0 = 3*ln(tau)
    public static IdealDerivatives Derivatives(double delta, double tau)
    {
        return new IdealDerivatives(3 * Math.Log(tau), 3 / tau, -3 / (tau * tau));
    }
}
