using System;

namespace EngineeringFluids.Helmholtz;
public static class IdealHelmholtzLogTau
{
    // alpha0 = 3*ln(tau)
    public static IdealDerivatives Derivatives(double delta, double tau) => Derivatives(delta, tau, Math.Log(tau));

    // logTau = Math.Log(tau), shared by every Helmholtz term of one state
    public static IdealDerivatives Derivatives(double delta, double tau, double logTau)
    {
        return new IdealDerivatives(3 * logTau, 3 / tau, -3 / (tau * tau));
    }
}
