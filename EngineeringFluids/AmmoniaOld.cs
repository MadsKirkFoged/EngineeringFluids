using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace EngineeringFluids;
public static class AmmoniaOld
{

    public static float alpha0(float delta, float tau)
    {

        return IdealGasHelmholtzLead(delta, tau) +
               IdealHelmholtzLogTau(delta, tau) +
               IdealHelmholtzPlanckEinstein(delta, tau);

    }

    public static float alphaR(float delta, float tau)
    { 
        return ResidualHelmholtzPower(delta, tau) +
               ResidualHelmholtzGaussian(delta, tau) +
               ResidualHelmholtzGaoB(delta, tau);
    
    }

    public static float alpha0_dTau(float delta, float tau)
    {
        var test1 = IdealGasHelmholtzLead_dTau(delta, tau);
        var test2 = IdealHelmholtzLogTau_dTau(delta, tau);
        var test3 = IdealHelmholtzPlanckEinstein_dTau(delta, tau);


        return IdealGasHelmholtzLead_dTau(delta, tau) +
           IdealHelmholtzLogTau_dTau(delta, tau) +
           IdealHelmholtzPlanckEinstein_dTau(delta, tau);

    }

    public static float alphaR_dDelta(float delta, float tau)
    {




        return ResidualHelmholtzPower_dDelta(delta, tau) +
               ResidualHelmholtzGaussian_dDelta(delta, tau) +
               ResidualHelmholtzGaoB_dDelta(delta, tau);

    }


    //TODO: Does not give exect the same as coolprop
    public static float alphaR_dTau(float delta, float tau)
    {
        var test1 = ResidualHelmholtzPower_dTau(delta, tau);
        //var test2 = ResidualHelmholtzPower_dTau2(delta, tau);
        //var test3 = ResidualHelmholtzPower_dTau3(delta, tau);
        //var test4 = ResidualHelmholtzPower_dTau4(delta, tau);

        var test4 = ResidualHelmholtzGaussian_dTau(delta, tau);
        var test5 = ResidualHelmholtzGaussian_dTau2(delta, tau);
        var test6 = ResidualHelmholtzGaussian_dTau3(delta, tau);
        var test7 = ResidualHelmholtzGaussian_dTau4(delta, tau);


        var test8 = ResidualHelmholtzGaoB_dTau(delta, tau);
        var test9 = ResidualHelmholtzGaoB_dTau2(delta, tau);
        var test10 = ResidualHelmholtzGaoB_dTau3(delta, tau);
        var test11 = ResidualHelmholtzGaoB_dTau4(delta, tau);

        return ResidualHelmholtzPower_dTau(delta, tau) +
               ResidualHelmholtzGaussian_dTau(delta, tau) +
               ResidualHelmholtzGaoB_dTau3(delta, tau);

    }


    //private static double IdealGasHelmholtzLead(double delta, double tau)
    //{
    //    double alpha0_IdealHelmholtzLead = Math.Log(delta) + -6.59406093943886 + 5.60101151987913 * tau;
    //    return alpha0_IdealHelmholtzLead;
    //}

    private static float IdealGasHelmholtzLead(float delta, float tau)
    {
        float alpha0_IdealHelmholtzLead = MathF.Log(delta) + -6.59406093943886f + 5.60101151987913f * tau;
        return alpha0_IdealHelmholtzLead;
    }

    private static float IdealGasHelmholtzLead_dTau(float delta, float tau)
    {
        // The derivative of the function with respect to tau
        float dalpha0_IdealHelmholtzLead_dTau = 5.60101151987913f;

        return dalpha0_IdealHelmholtzLead_dTau;
    }

    //private static double IdealHelmholtzLogTau(double delta, double tau)
    //{
    //    double alpha0 = 3 * Math.Log(tau);
    //    return alpha0;
    //}

    private static float IdealHelmholtzLogTau(float delta, float tau)
    {
        float alpha0 = 3 * MathF.Log(tau);
        return alpha0;
    }

    private static float IdealHelmholtzLogTau_dTau(float delta, float tau)
    {
        // The derivative of the function with respect to tau
        float dalpha0_dTau = 3 / tau;

        return dalpha0_dTau;
    }

    private static readonly double[] n_PE = [2.224, 3.148, 0.9579 ];
    private static readonly double[] t_PE = [-4.0585856593352405, -9.776605187888352, -17.829667620080876];

    private static readonly double[] c_PE = [1, 1, 1];
    private static readonly double[] d_PE = [-1, -1, -1];

    private static readonly float[] n_PEF = [2.224f, 3.148f, 0.9579f];
    private static readonly float[] t_PEF = [-4.0585856593352405f, -9.776605187888352f, -17.829667620080876f];

    private static readonly float[] c_PEF = [1, 1, 1];
    private static readonly float[] d_PEF = [-1, -1, -1];



    //private static double IdealHelmholtzPlanckEinstein(double delta, double tau)
    //{
    //    double alpha0_PE = 0.0;
    //    for (int i = 0; i < n_PE.Length; i++)
    //    {
    //        alpha0_PE += n_PE[i] * Math.Log(c_PE[i] + d_PE[i] * Math.Exp(t_PE[i] * tau));
    //        //alpha0_PE += n_PE[i] * Math.Log(-1 * Math.Exp(t_PE[i] * tau));


    //    }

    //    return alpha0_PE;
    //}

    public static float IdealHelmholtzPlanckEinstein(float delta, float tau)
    {
        float alpha0_PE = 0.0f;
        for (int i = 0; i < n_PEF.Length; i++)
        {
            alpha0_PE += n_PEF[i] * MathF.Log(c_PEF[i] + d_PEF[i] * MathF.Exp(t_PEF[i] * tau));
        }


        return alpha0_PE;
    }

    public static float IdealHelmholtzPlanckEinstein_dTau(float delta, float tau)
    {
        float dalpha0_PE_dTau = 0.0f;
        for (int i = 0; i < n_PEF.Length; i++)
        {
            float numerator = n_PEF[i] * d_PEF[i] * t_PEF[i] * MathF.Exp(t_PEF[i] * tau);
            float denominator = c_PEF[i] + d_PEF[i] * MathF.Exp(t_PEF[i] * tau);
            dalpha0_PE_dTau += numerator / denominator;
        }

        return dalpha0_PE_dTau;
    }   

    private static readonly float[] n_PESpan = [2.224f, 3.148f, 0.9579f,0f, 0f, 0f, 0f, 0f];
    private static readonly float[] t_PESpan = [-4.0585856593352405f, -9.776605187888352f, -17.829667620080876f, 0f, 0f, 0f, 0f, 0f];

    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
    //public static float IdealHelmholtzPlanckEinsteinCIMD(float delta, float tau)
    //{
    //    Span<float> test = stackalloc float[8];
    //    test.Fill(tau);

    //    FastMath.IdealHelmholtzPlanckEinstein(test, n_PESpan, t_PESpan);

    //    return test[0];
    //}






    private static readonly float[] n_RF = { 0.006132232f, 1.7395866f, -2.2261792f, -0.30127553f, 0.08967023f, -0.076387037f, -0.84063963f, -0.27026327f };
    private static readonly double[] n_R = { 0.006132232, 1.7395866, -2.2261792, -0.30127553, 0.08967023, -0.076387037, -0.84063963, -0.27026327 };
    private static readonly int[] d_R = { 4, 1, 1, 2, 3, 3, 2, 3 };
    private static readonly float[] d_RF = { 4, 1, 1, 2, 3, 3, 2, 3 };
    private static readonly double[] t_R = { 1.0, 0.382, 1.0, 1.0, 0.677, 2.915, 3.51, 1.063 };
    private static readonly float[] t_RF = { 1.0f, 0.382f, 1.0f, 1.0f, 0.677f, 2.915f, 3.51f, 1.063f };
    private static readonly int[] I_R = { 0, 0, 0, 0, 0, 2, 2, 1 };
    private static readonly float[] I_Rf = { 0f, 0f, 0f, 0f, 0f, 2f, 2f, 1f };

    private static double Pow(double x, int y)
    {
        var result = 1.0;

        for (int i = 0; i < y; i++)
        {
            result *= x;
        }
        return result;
    }

    private static float Pow(float x, int y)
    {
        var result = 1f;

        for (int i = 0; i < y; i++)
        {
            result *= x;
        }
        return result;
    }

    private static double Pow(double x, double y)
    {
        int i = (int)x;

        if (x == i)
        {
            return Pow(x, i);
        }
        else
        {
            return Math.Pow(x, y);
        }
    }

    private static float Pow(float x, float y)
    {
        int i = (int)x;

        if (x == i)
        {
            return Pow(x, i);
        }
        else
        {
            return MathF.Pow(x, y);
        }
    }

    public static double[] Multiply(double[] left, double[] right, double[] extra)
    {
        if (left == null)
            throw new ArgumentNullException(nameof(left));
        if (right == null)
            throw new ArgumentNullException(nameof(right));
        if (extra == null)
            throw new ArgumentNullException(nameof(extra));
        if (left.Length != right.Length || right.Length != extra.Length)
            throw new ArgumentException("All arrays must have the same length");

        int length = left.Length;
        double[] result = new double[length];

        int remaining = length % Vector<double>.Count;

        for (int i = 0; i < length - remaining; i += Vector<double>.Count)
        {
            var v1 = new Vector<double>(left, i);
            var v2 = new Vector<double>(right, i);
            var v3 = new Vector<double>(extra, i);
            (v1 * v2 * v3).CopyTo(result, i);
        }

        for (int i = length - remaining; i < length; i++)
        {
            result[i] = left[i] * right[i] * extra[i];
        }

        return result;
    }


    private static float[] PowVec(float x, float[] right)
    {

        int length = right.Length;
        float[] result = new float[length];
        float[] left = new float[length];


        // Get the number of elements that can't be processed in the vector
        // NOTE: Vector<T>.Count is a JIT time constant and will get optimized accordingly
        int remaining = length % Vector<float>.Count;

        for (int i = 0; i < length - remaining; i += Vector<float>.Count)
        {
            var v1 = new Vector<float>(left, i);
            var v2 = new Vector<float>(right, i);
            (v1 * v2).CopyTo(result, i);
        }

        for (int i = length - remaining; i < length; i++)
        {
            result[i] = left[i] + right[i];
        }

        return result;
    }

    private static float[] MultiplyVec(float[] first, float[] second, float[] third)
    {

        int length = first.Length;
        float[] result = new float[length];


        // Get the number of elements that can't be processed in the vector
        // NOTE: Vector<T>.Count is a JIT time constant and will get optimized accordingly
        int remaining = length % Vector<float>.Count;

        for (int i = 0; i < length - remaining; i += Vector<float>.Count)
        {
            var v1 = new Vector<float>(first, i);
            var v2 = new Vector<float>(second, i);
            var v3 = new Vector<float>(third, i);
            (v1 * v2 * v3).CopyTo(result, i);

        }

        for (int i = length - remaining; i < length; i++)
        {
            result[i] = first[i] * second[i] * third[i];
        }

        return result;
    }


  //  ResidualHelmholtzPower:
  //α^r = { 
  //      ∑_i n_i δ^{d_i}τ^{t_i}                    if l_i = 0
  //      ∑_i n_i δ^{d_i} τ^{ t_i} exp(-δ^{ l_i})   if l_i ≠ 0
  //      }

    public static float ResidualHelmholtzPower(float delta, float tau)
    {
        float alphaRP = 0.0f;


        for (int i = 0; i < n_R.Length; i++)
        {
            float term = n_RF[i] * Pow(delta, d_R[i]) * Pow(tau, t_RF[i]);

            if (I_R[i] != 0)
            {
                term *= MathF.Exp(-Pow(delta, I_R[i]));
            }
            alphaRP += term;
        }

        return alphaRP;
    }

    
    public static float ResidualHelmholtzPower_dDelta(float delta, float tau)
    {
        float dalphaRP_dDelta = 0.0f;

        for (int i = 0; i < n_RF.Length; i++)
        {
            // Base term without the exponential component
            float baseTerm = n_RF[i] * Pow(tau, t_RF[i]);
            // Derivative of the power term
            float powerTerm = d_R[i] * Pow(delta, d_R[i] - 1);
            // Initial term calculation
            float term = baseTerm * powerTerm;

            if (I_R[i] != 0)
            {
                // Original term with exponential: n_RF[i] * delta^(d_R[i]) * tau^(t_RF[i]) * exp(-delta^(I_R[i]))
                // Derivative: 
                // d_R[i] * n_RF[i] * delta^(d_R[i] - 1) * tau^(t_RF[i]) * exp(-delta^(I_R[i]))
                // - I_R[i] * n_RF[i] * delta^(d_R[i] + I_R[i] - 1) * tau^(t_RF[i]) * exp(-delta^(I_R[i]))
                float expTerm = MathF.Exp(-Pow(delta, I_R[i]));
                float expDerivativeTerm = -I_R[i] * Pow(delta, I_R[i] - 1) * expTerm;

                term = baseTerm * (powerTerm * expTerm + Pow(delta, d_R[i]) * expDerivativeTerm);
            }

            dalphaRP_dDelta += term;
        }

        return dalphaRP_dDelta;
    }

    public static float ResidualHelmholtzPower_dTau(float delta, float tau)
    {
        float dalphaRP_dTau = 0.0f;

        for (int i = 0; i < n_R.Length; i++)
        {
            float term = n_RF[i] * MathF.Pow(delta, d_R[i]) * MathF.Pow(tau, t_RF[i]);

            if (I_R[i] != 0)
            {
                term *= MathF.Exp(-MathF.Pow(delta, I_R[i]));
            }

            // Differentiate term with respect to tau
            float dterm_dTau = term * t_RF[i] / tau;

            dalphaRP_dTau += dterm_dTau;
        }

        return dalphaRP_dTau;
    }    


    private static readonly float[] n_Rg = { 6.212578f, -5.7844357f, 2.4817542f, -2.3739168f, 0.01493697f, -3.7749264f, 0.0006254348f, -1.7359e-05f, -0.13462033f, 0.07749072839f };
   private static readonly float[] d_Rg = { 1, 1, 1, 2, 2, 1, 3, 3, 1, 1 };
   private static readonly float[] t_Rg = { 0.655f, 1.3f, 3.1f, 1.4395f, 1.623f, 0.643f, 1.13f, 4.5f, 1.0f, 4.0f };
   private static readonly float[] eta_Rg = { 0.42776f, 0.6424f, 0.8175f, 0.7995f, 0.91f, 0.3574f, 1.21f, 4.14f, 22.56f, 22.68f };
   private static readonly float[] epsilon_Rg = { -0.0726f, -0.1274f, 0.7527f, 0.57f, 2.2f, -0.243f, 2.96f, 3.02f, 0.9574f, 0.9576f };
   private static readonly float[] beta_Rg = { 1.708f, 1.4865f, 2.0915f, 2.43f, 0.488f, 1.1f, 0.85f, 1.14f, 945.64f, 993.85f };
   private static readonly float[] gamma_Rg = { 1.036f, 1.2777f, 1.083f, 1.2906f, 0.928f, 0.934f, 0.919f, 1.852f, 1.05897f, 1.05277f };

    public static float ResidualHelmholtzGaussian(float delta, float tau)
    {
        // Calculate alphaR
        float alphaRg = 0.0f;
        for (int i = 0; i < n_R.Length; i++)
        {
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);
            float tauTerm = MathF.Pow(tau, t_Rg[i]);
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            float term = n_Rg[i] * deltaTerm * tauTerm * expTerm;
            alphaRg += term;
        }
        return alphaRg;
    }

    public static float ResidualHelmholtzGaussian_dDelta(float delta, float tau)
    {
        float dalphaRg_dDelta = 0.0f;

        for (int i = 0; i < n_R.Length; i++)
        {
            // Original deltaTerm: delta^(d_Rg[i])
            // Derivative of deltaTerm: d_Rg[i] * delta^(d_Rg[i] - 1)
            float deltaTerm = d_Rg[i] * MathF.Pow(delta, d_Rg[i] - 1);

            // Original tauTerm: tau^(t_Rg[i])
            float tauTerm = MathF.Pow(tau, t_Rg[i]);

            // Original expTerm: exp(-eta_Rg[i] * (delta - epsilon_Rg[i])^2 - beta_Rg[i] * (tau - gamma_Rg[i])^2)
            // Derivative of expTerm with respect to delta:
            // -2 * eta_Rg[i] * (delta - epsilon_Rg[i]) * exp(-eta_Rg[i] * (delta - epsilon_Rg[i])^2 - beta_Rg[i] * (tau - gamma_Rg[i])^2)
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            float expTermDerivative = -2 * eta_Rg[i] * (delta - epsilon_Rg[i]) * expTerm;

            // Combined term:
            // term = n_Rg[i] * (deltaTerm * tauTerm * expTermDerivative + d_Rg[i] * delta^(d_Rg[i] - 1) * tauTerm * expTerm)
            float term = n_Rg[i] * (deltaTerm * tauTerm * expTerm + MathF.Pow(delta, d_Rg[i]) * tauTerm * expTermDerivative);

            dalphaRg_dDelta += term;
        }

        return dalphaRg_dDelta;
    }

    public static float ResidualHelmholtzGaussian_dTau(float delta, float tau)
    {
        float dalphaRg_dTau = 0.0f;
        for (int i = 0; i < n_Rg.Length; i++)
        {
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);
            float tauTerm = MathF.Pow(tau, t_Rg[i]);
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            float term = n_Rg[i] * deltaTerm * tauTerm * expTerm;

            // Differentiate term with respect to tau
            float dterm_dTau = term * (t_Rg[i] / tau - 2 * beta_Rg[i] * (tau - gamma_Rg[i]));

            dalphaRg_dTau += dterm_dTau;
        }
        return dalphaRg_dTau;
    }

    public static float ResidualHelmholtzGaussian_dTau2(float delta, float tau)
    {
        float dalphaRg_dTau = 0.0f;
        for (int i = 0; i < n_Rg.Length; i++)
        {
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);
            float tauTerm = MathF.Pow(tau, t_Rg[i]);
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            float term = n_Rg[i] * deltaTerm * tauTerm * expTerm;

            // Differentiate term with respect to tau
            float dterm_dTau = term * ((t_Rg[i] / tau) - (2 * beta_Rg[i] * (tau - gamma_Rg[i])));

            dalphaRg_dTau += dterm_dTau;
        }
        return dalphaRg_dTau;
    }

    public static float ResidualHelmholtzGaussian_dTau3(float delta, float tau)
    {
        float dalphaRg_dTau = 0.0f;

        for (int i = 0; i < n_R.Length; i++)
        {
            // Original deltaTerm: delta^(d_Rg[i])
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);

            // Derivative of tauTerm: t_Rg[i] * tau^(t_Rg[i] - 1)
            float tauTerm = t_Rg[i] * MathF.Pow(tau, t_Rg[i] - 1);

            // Original expTerm: exp(-eta_Rg[i] * (delta - epsilon_Rg[i])^2 - beta_Rg[i] * (tau - gamma_Rg[i])^2)
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            // Derivative of expTerm with respect to tau:
            float expTermDerivative = -2 * beta_Rg[i] * (tau - gamma_Rg[i]) * expTerm;

            // Combined term:
            // term = n_Rg[i] * (deltaTerm * tauTerm * expTerm + deltaTerm * tau^(t_Rg[i]) * expTermDerivative)
            float term = n_Rg[i] * (deltaTerm * tauTerm * expTerm + deltaTerm * MathF.Pow(tau, t_Rg[i]) * expTermDerivative);

            dalphaRg_dTau += term;
        }

        return dalphaRg_dTau;
    }


    public static float ResidualHelmholtzGaussian_dTau4(float delta, float tau)
    {
        float dalphaRg_dTau = 0.0f;

        for (int i = 0; i < n_R.Length; i++)
        {
            // deltaTerm: delta^(d_Rg[i])
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);

            // Derivative of tauTerm: t_Rg[i] * tau^(t_Rg[i] - 1)
            float tauTermDerivative = t_Rg[i] * MathF.Pow(tau, t_Rg[i] - 1);

            // Original tauTerm: tau^(t_Rg[i])
            float tauTerm = MathF.Pow(tau, t_Rg[i]);

            // Original expTerm: exp(-eta_Rg[i] * (delta - epsilon_Rg[i])^2 - beta_Rg[i] * (tau - gamma_Rg[i])^2)
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));

            // Derivative of expTerm with respect to tau:
            float expTermDerivative = -2 * beta_Rg[i] * (tau - gamma_Rg[i]) * expTerm;

            // Combined term:
            // term = n_Rg[i] * deltaTerm * (tauTermDerivative * expTerm + tauTerm * expTermDerivative)
            float term = n_Rg[i] * deltaTerm * (tauTermDerivative * expTerm + tauTerm * expTermDerivative);

            dalphaRg_dTau += term;
        }

        return dalphaRg_dTau;
    }


    // Coefficients for Residual Helmholtz GaoB (alphaR)
    private static readonly float[] b = { 1.244f, 0.6826f };
   private static readonly float[] beta = { 0.3696f, 0.2962f };
   private static readonly float[] d = { 1, 1 };
   private static readonly float[] epsilon = { 0.4478f, 0.44689f };
   private static readonly float[] eta = { -2.8452f, -2.8342f };
   private static readonly float[] gamma = { 1.108f, 1.313f };
   private static readonly float[] n = { -1.6909858f, 0.93739074f };
   private static readonly float[] t = { 4.3315f, 4.015f };

    public static float ResidualHelmholtzGaoB(float delta, float tau)
    {
        // Calculate alphaR
        float alphaRb = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float Fdelta = MathF.Pow(delta, d[i]) * MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));
            float term = n[i] * Ftau * Fdelta;
            alphaRb += term;
        }
        return alphaRb;
    }

    public static float ResidualHelmholtzGaoB_dDelta(float delta, float tau)
    {
        float dalphaRb_dDelta = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {
            // Original Ftau term: tau^(t[i]) * exp(1 / (b[i] + beta[i] * (gamma[i] - tau)^2))
            // Ftau does not depend on delta, so its derivative with respect to delta is 0
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));

            // Original Fdelta term: delta^(d[i]) * exp(eta[i] * (delta - epsilon[i])^2)
            // Derivative of delta^(d[i]) with respect to delta: d[i] * delta^(d[i] - 1)
            // Derivative of exp(eta[i] * (delta - epsilon[i])^2) with respect to delta:
            // 2 * eta[i] * (delta - epsilon[i]) * exp(eta[i] * (delta - epsilon[i])^2)
            float deltaTerm = d[i] * MathF.Pow(delta, d[i] - 1);
            float expTerm = MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));
            float expTermDerivative = 2 * eta[i] * (delta - epsilon[i]) * expTerm;

            // Combined term for derivative:
            // term = n[i] * (Ftau * (deltaTerm * expTerm + delta^(d[i]) * expTermDerivative))
            float term = n[i] * Ftau * (deltaTerm * expTerm + MathF.Pow(delta, d[i]) * expTermDerivative);

            dalphaRb_dDelta += term;
        }

        return dalphaRb_dDelta;
    }


    public static float ResidualHelmholtzGaoB_dTau(float delta, float tau)
    {
        float dalphaRb_dTau = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float Fdelta = MathF.Pow(delta, d[i]) * MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));

            // Calculate the derivative of Ftau with respect to tau
            float dFtau_dTau = MathF.Pow(tau, t[i] - 1) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2))) *
                               (t[i] - 2 * beta[i] * (gamma[i] - tau) / MathF.Pow(b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2), 2));

            float term = n[i] * Fdelta * dFtau_dTau;

            dalphaRb_dTau += term;
        }
        return dalphaRb_dTau;
    }

    public static float ResidualHelmholtzGaoB_dTau2(float delta, float tau)
    {
        float dalphaRb_dTau = 0.0f;
        for (int i = 0; i < n.Length; i++)
        {
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float Fdelta = MathF.Pow(delta, d[i]) * MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));

            // Calculate the derivative of Ftau with respect to tau
            float tauTerm = MathF.Pow(tau, t[i] - 1);
            float expTerm = MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float dFtau_dTau = tauTerm * expTerm *
                               (t[i] - 2 * beta[i] * (gamma[i] - tau) / MathF.Pow(b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2), 2));

            float term = n[i] * Fdelta * dFtau_dTau;

            dalphaRb_dTau += term;
        }
        return dalphaRb_dTau;
    }

    public static float ResidualHelmholtzGaoB_dTau3(float delta, float tau)
    {
        float dalphaRb_dTau = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {
            // Original Ftau term: tau^(t[i]) * exp(1 / (b[i] + beta[i] * (gamma[i] - tau)^2))
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));

            // Derivative of Ftau term with respect to tau
            float expTerm = MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float expTermDerivative = expTerm * (-2 * beta[i] * (gamma[i] - tau) / MathF.Pow(b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2), 2));

            float FtauDerivative = t[i] * MathF.Pow(tau, t[i] - 1) * expTerm + MathF.Pow(tau, t[i]) * expTermDerivative;

            // Original Fdelta term: delta^(d[i]) * exp(eta[i] * (delta - epsilon[i])^2)
            float Fdelta = MathF.Pow(delta, d[i]) * MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));

            // Combined term for derivative
            float term = n[i] * FtauDerivative * Fdelta;

            dalphaRb_dTau += term;
        }

        return dalphaRb_dTau;
    }

    public static float ResidualHelmholtzGaoB_dTau4(float delta, float tau)
    {
        float dalphaRb_dTau = 0.0f;

        for (int i = 0; i < n.Length; i++)
        {
            // Original Ftau term: tau^(t[i]) * exp(1 / (b[i] + beta[i] * (gamma[i] - tau)^2))
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));

            // Derivative of Ftau term with respect to tau
            float expTerm = MathF.Exp(1.0f / (b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2)));
            float expTermDerivative = -2.0f * beta[i] * (gamma[i] - tau) / MathF.Pow(b[i] + beta[i] * MathF.Pow(gamma[i] - tau, 2), 2);

            float FtauDerivative = t[i] * MathF.Pow(tau, t[i] - 1) * expTerm + MathF.Pow(tau, t[i]) * expTermDerivative * expTerm;

            // Original Fdelta term: delta^(d[i]) * exp(eta[i] * (delta - epsilon[i])^2)
            float Fdelta = MathF.Pow(delta, d[i]) * MathF.Exp(eta[i] * MathF.Pow(delta - epsilon[i], 2));

            // Combined term for derivative
            float term = n[i] * FtauDerivative * Fdelta;

            dalphaRb_dTau += term;
        }

        return dalphaRb_dTau;
    }





}
