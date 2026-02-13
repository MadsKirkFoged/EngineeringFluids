using EngineeringMath.VectorMath.Add;
using EngineeringMath.VectorMath.Exp;
using EngineeringMath.VectorMath.Pow;
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
public static class AmmoniaFast
{
    public static float alpha(float delta, float tau)
    {

        return alpha0(delta, tau) + alphaR(delta, tau);

    }



    //[MethodImpl(MethodImplOptions.AggressiveInlining)]
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


    private static float IdealGasHelmholtzLead(float delta, float tau)
    {
        float alpha0_IdealHelmholtzLead = MathF.Log(delta) + -6.59406093943886f + 5.60101151987913f * tau;
        return alpha0_IdealHelmholtzLead;
    }

    private static float IdealHelmholtzLogTau(float delta, float tau)
    {
        float alpha0 = 3 * MathF.Log(tau);
        return alpha0;
    }

    private static readonly float[] n_PEF = [2.224f, 3.148f, 0.9579f];
    private static readonly float[] t_PEF = [-4.0585856593352405f, -9.776605187888352f, -17.829667620080876f];

    private static readonly float[] c_PEF = [1, 1, 1];
    private static readonly float[] d_PEF = [-1, -1, -1];


    public static float IdealHelmholtzPlanckEinstein(float delta, float tau)
    {
        float alpha0_PE = 0.0f;
        for (int i = 0; i < n_PEF.Length; i++)
        {
            alpha0_PE += n_PEF[i] * MathF.Log(c_PEF[i] + d_PEF[i] * MathF.Exp(t_PEF[i] * tau));
        }

        return alpha0_PE;
    }


    private static readonly float[] n_R = { 0.006132232f, 1.7395866f, -2.2261792f, -0.30127553f, 0.08967023f, -0.076387037f, -0.84063963f, -0.27026327f };  
    private static readonly int[] d_R = { 4, 1, 1, 2, 3, 3, 2, 3 };
    private static readonly float[] t_R = { 1.0f, 0.382f, 1.0f, 1.0f, 0.677f, 2.915f, 3.51f, 1.063f };
    private static readonly int[] I_R = { 0, 0, 0, 0, 0, 2, 2, 1 };


    //    | Method  | Mean      | Error    | StdDev   | Gen0   | Allocated |
    //    |-------- |----------:|---------:|---------:|-------:|----------:|
    //    | org     | 123.51 ns | 1.317 ns | 1.028 ns |      - |         - |
    //    | alpha02 | 67.60 ns  | 0.946 ns | 0.790 ns |      - |         - |



    //Is the fastest
    public static float ResidualHelmholtzPower(float delta, float tau)
    {
        Span<float> deltaArray = stackalloc float[8];
        deltaArray.Fill(delta);

        Span<float> tauArray = stackalloc float[8];
        tauArray.Fill(tau);

        Power8.Pow(deltaArray, d_R);
        Power8.Pow(tauArray, t_R);

        Multiply8.MultiplyFixed8Reuse(deltaArray, tauArray, n_R);


        for (int i = 0; i < deltaArray.Length; i++)
        {
            if (I_R[i] != 0)
            {
                deltaArray[i] *= MathF.Exp(-Power.Pow(delta, I_R[i]));
            }
        }


        Sum8.Sum(deltaArray);


        return deltaArray[0];
    }





    private static readonly float[] n_Rg = { 6.212578f, -5.7844357f, 2.4817542f, -2.3739168f, 0.01493697f, -3.7749264f, 0.0006254348f, -1.7359e-05f, -0.13462033f, 0.07749072839f };
    private static readonly int[] d_Rg = { 1, 1, 1, 2, 2, 1, 3, 3, 1, 1 };
    private static readonly float[] t_Rg = { 0.655f, 1.3f, 3.1f, 1.4395f, 1.623f, 0.643f, 1.13f, 4.5f, 1.0f, 4.0f };
    private static readonly float[] eta_Rg = { -0.42776f, -0.6424f, -0.8175f, -0.7995f, -0.91f, -0.3574f, -1.21f, -4.14f, -22.56f, -22.68f };
    private static readonly float[] epsilon_Rg = { -0.0726f, -0.1274f, 0.7527f, 0.57f, 2.2f, -0.243f, 2.96f, 3.02f, 0.9574f, 0.9576f };
    private static readonly float[] beta_Rg = { 1.708f, 1.4865f, 2.0915f, 2.43f, 0.488f, 1.1f, 0.85f, 1.14f, 945.64f, 993.85f };
    private static readonly float[] gamma_Rg = { 1.036f, 1.2777f, 1.083f, 1.2906f, 0.928f, 0.934f, 0.919f, 1.852f, 1.05897f, 1.05277f };

    public static double ResidualHelmholtzGaussian(double delta, double tau)
    {
        // Calculate alphaR
        double alphaRg = 0.0;
        for (int i = 0; i < n_R.Length; i++)
        {
            double deltaTerm = Math.Pow(delta, d_Rg[i]);
            double tauTerm = Math.Pow(tau, t_Rg[i]);
            double expTerm = Math.Exp(-eta_Rg[i] * Math.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * Math.Pow(tau - gamma_Rg[i], 2));
            double term = n_Rg[i] * deltaTerm * tauTerm * expTerm;
            alphaRg += term;
        }
        return alphaRg;
    }


    public static float ResidualHelmholtzGaussian(float delta, float tau)
    {

        Span<float> deltaArray = stackalloc float[8];
        deltaArray.Fill(delta);

        Span<float> deltaArray2 = stackalloc float[8];
        deltaArray2.Fill(delta);

        Span<float> tauArray = stackalloc float[8];
        tauArray.Fill(tau);

        Span<float> tauArray2 = stackalloc float[8];
        tauArray2.Fill(tau);

        Power8.Pow(deltaArray, d_Rg);
        Power8.Pow(tauArray, t_Rg);

        //double expTerm = Math.Exp(-eta_Rg[i] * Math.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i];
        Subtraction8.SubFixed8Reuse(tauArray2, gamma_Rg);
        Power8.Pow(tauArray2, 2);
        Multiply8.MultiplyFixed8Reuse(tauArray2, beta_Rg);

        Subtraction8.SubFixed8Reuse(deltaArray2, epsilon_Rg);
        Power8.Pow(deltaArray2, 2);
        Multiply8.MultiplyFixed8Reuse(deltaArray2, eta_Rg);

        Subtraction8.SubFixed8Reuse(deltaArray2, tauArray2);
        ExpClass.Exp8(deltaArray2);

        Multiply8.MultiplyFixed8Reuse(deltaArray, tauArray, deltaArray2, n_Rg);

        Sum8.Sum(deltaArray);



        float alphaRg = 0.0f;
        for (int i = 8; i < n_R.Length; i++)
        {
            float deltaTerm = MathF.Pow(delta, d_Rg[i]);
            float tauTerm = MathF.Pow(tau, t_Rg[i]);
            float expTerm = MathF.Exp(-eta_Rg[i] * MathF.Pow(delta - epsilon_Rg[i], 2) - beta_Rg[i] * MathF.Pow(tau - gamma_Rg[i], 2));
            float term = n_Rg[i] * deltaTerm * tauTerm * expTerm;
            alphaRg += term;
        }


        return alphaRg + deltaArray[0];
    }





    private static readonly float[] b = { 1.244f, 0.6826f };
    private static readonly float[] beta = { 0.3696f, 0.2962f };
    private static readonly int[] d = { 1, 1 };
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
            float Ftau = MathF.Pow(tau, t[i]) * MathF.Exp(1.0f / (b[i] + beta[i] * Power.Pow(gamma[i] - tau, 2)));
            float Fdelta = Power.Pow(delta, d[i]) * MathF.Exp(eta[i] * Power.Pow(delta - epsilon[i], 2));
            float term = n[i] * Ftau * Fdelta;
            alphaRb += term;
        }
        return alphaRb;
    }


}
