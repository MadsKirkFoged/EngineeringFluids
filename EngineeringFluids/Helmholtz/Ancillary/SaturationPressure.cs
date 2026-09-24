using System;
using System.Runtime.CompilerServices;

/// <summary>
/// Fast pressure approximation generated from TemperatureThenPressure.txt.
/// Two-part piecewise quadratic:
///   - T <= 403°C: 2°C segments, quadratic in normalized u
///   - T >  403°C: 0.5°C segments, quadratic in normalized u
///
/// Form per segment: P(T) ~= c0 + c1*u + c2*u^2, with u in [-1,1]
///
/// Segment mapping is O(1) (no binary search), suited for hot paths.
/// </summary>
public static class SaturationPressureFast
{
    // Domain (from your old system)
    public const float Tmin = 195.495f;
    public const float Tmax = 405.56f;
    public const float Tsplit = 403.0f;

    // ===== Low region (<= 403): 2°C segments, half-width = 1 =====
    private const float SegW_Lo = 2.0f;
    private const float InvSegW_Lo = 0.5f;         // 1/2
    private const float Half_Lo = 1.0f;            // SegW/2
    private const float TAnchor_Lo = 194.995000000f;
    private const int NSeg_Lo = 104;

    // ===== High region (> 403): 0.5°C segments, half-width = 0.25 =====
    private const float SegW_Hi = 0.5f;
    private const float InvSegW_Hi = 2.0f;         // 1/0.5
    private const float Half_Hi = 0.25f;           // SegW/2
    private const float TAnchor_Hi = 402.555000000f;
    private const int NSeg_Hi = 6;

    /// <summary>
    /// Returns pressure (Pa) for a given temperature (°C).
    /// </summary>
    //[MethodImpl(MethodImplOptions.AggressiveInlining | MethodImplOptions.AggressiveOptimization)]
    public static float Pressure(float t)
    {
        // Clamp to domain
        if (t < Tmin)
            t = Tmin;
        else if (t > Tmax)
            t = Tmax;

        if (t <= Tsplit)
        {
            // k = floor((t - TAnchor_Lo) / 2)
            int k = (int)((t - TAnchor_Lo) * InvSegW_Lo);
            if ((uint)k >= (uint)NSeg_Lo)
                k = (k < 0) ? 0 : (NSeg_Lo - 1);

            // mid = TAnchor_Lo + (k + 0.5) * 2
            float mid = TAnchor_Lo + (k + 0.5f) * SegW_Lo;
            float u = (t - mid) / Half_Lo; // Half_Lo=1 => u = t - mid

            int o = k * 3;
            float c0 = CoeffLow[o + 0];
            float c1 = CoeffLow[o + 1];
            float c2 = CoeffLow[o + 2];

            // Horner: (c2*u + c1)*u + c0
            //If you want FMA: 
               float y = MathF.FusedMultiplyAdd(c2, u, c1);
            return MathF.FusedMultiplyAdd(y, u, c0);  // see docs [1](https://learn.microsoft.com/en-us/dotnet/api/system.math.fusedmultiplyadd?view=net-10.0)
            //return (c2 * u + c1) * u + c0;
        }
        else
        {
            // k = floor((t - TAnchor_Hi) / 0.5)
            int k = (int)((t - TAnchor_Hi) * InvSegW_Hi);
            if ((uint)k >= (uint)NSeg_Hi)
                k = (k < 0) ? 0 : (NSeg_Hi - 1);

            // mid = TAnchor_Hi + (k + 0.5) * 0.5
            float mid = TAnchor_Hi + (k + 0.5f) * SegW_Hi;
            float u = (t - mid) / Half_Hi; // Half_Hi=0.25 -> u = 4*(t-mid)

            int o = k * 3;
            float c0 = CoeffHigh[o + 0];
            float c1 = CoeffHigh[o + 1];
            float c2 = CoeffHigh[o + 2];

            float y = MathF.FusedMultiplyAdd(c2, u, c1);
            return MathF.FusedMultiplyAdd(y, u, c0);

            //return (c2 * u + c1) * u + c0;
        }
    }

    // Coeff layout: [seg0 c0,c1,c2, seg1 c0,c1,c2, ...]
    private static ReadOnlySpan<float> CoeffLow => new float[]
    {
        6.30270459e+3f, 5.02345642e+2f, 1.73336773e+1f, 7.37828906e+3f, 5.74783325e+2f, 1.90657024e+1f,
        8.60650195e+3f, 6.55273376e+2f, 2.11934452e+1f, 1.00043926e+4f, 7.44605225e+2f, 2.34866619e+1f,
        1.15903145e+4f, 8.43452454e+2f, 2.59511623e+1f, 1.33839883e+4f, 9.52511353e+2f, 2.85924454e+1f,
        1.54065527e+4f, 1.07249963e+3f, 3.14156685e+1f, 1.76805977e+4f, 1.20415479e+3f, 3.44256325e+1f,
        2.02302109e+4f, 1.34823267e+3f, 3.76267548e+1f, 2.30810059e+4f, 1.50550623e+3f, 4.10230675e+1f,
        2.62601602e+4f, 1.67676355e+3f, 4.46181831e+1f, 2.97964395e+4f, 1.86280652e+3f, 4.84153175e+1f,
        3.37202266e+4f, 2.06444873e+3f, 5.24172592e+1f, 3.80635430e+4f, 2.28251465e+3f, 5.66263847e+1f,
        4.28600625e+4f, 2.51783691e+3f, 6.10446510e+1f, 4.81451445e+4f, 2.77125537e+3f, 6.56736069e+1f,
        5.39558203e+4f, 3.04361523e+3f, 7.05143967e+1f, 6.03308203e+4f, 3.33576514e+3f, 7.55677872e+1f,
        6.73105781e+4f, 3.64855664e+3f, 8.08341446e+1f, 7.49372266e+4f, 3.98284180e+3f, 8.63134766e+1f,
        8.32546094e+4f, 4.33947217e+3f, 9.20054703e+1f, 9.23082578e+4f, 4.71929639e+3f, 9.79094543e+1f,
        1.02145414e+5f, 5.12316113e+3f, 1.04024460e+2f, 1.12815008e+5f, 5.55190820e+3f, 1.10349243e+2f,
        1.24367625e+5f, 6.00637354e+3f, 1.16882278e+2f, 1.36855547e+5f, 6.48738623e+3f, 1.23621803e+2f,
        1.50332688e+5f, 6.99576904e+3f, 1.30565842e+2f, 1.64854594e+5f, 7.53233545e+3f, 1.37712204e+2f,
        1.80478453e+5f, 8.09788965e+3f, 1.45058533e+2f, 1.97263031e+5f, 8.69322656e+3f, 1.52602325e+2f,
        2.15268688e+5f, 9.31913184e+3f, 1.60340927e+2f, 2.34557312e+5f, 9.97637793e+3f, 1.68271591e+2f,
        2.55192391e+5f, 1.06657285e+4f, 1.76391449e+2f, 2.77238844e+5f, 1.13879326e+4f, 1.84697586e+2f,
        3.00763156e+5f, 1.21437314e+4f, 1.93187057e+2f, 3.25833219e+5f, 1.29338516e+4f, 2.01856812e+2f,
        3.52518406e+5f, 1.37590078e+4f, 2.10703842e+2f, 3.80889500e+5f, 1.46199033e+4f, 2.19725128e+2f,
        4.11018688e+5f, 1.55172295e+4f, 2.28917633e+2f, 4.42979469e+5f, 1.64516641e+4f, 2.38278397e+2f,
        4.76846750e+5f, 1.74238750e+4f, 2.47804474e+2f, 5.12696750e+5f, 1.84345176e+4f, 2.57492981e+2f,
        5.50606938e+5f, 1.94842363e+4f, 2.67341125e+2f, 5.90656188e+5f, 2.05736641e+4f, 2.77346161e+2f,
        6.32924438e+5f, 2.17034219e+4f, 2.87505493e+2f, 6.77493062e+5f, 2.28741230e+4f, 2.97816589e+2f,
        7.24444500e+5f, 2.40863691e+4f, 3.08277100e+2f, 7.73862438e+5f, 2.53407559e+4f, 3.18884766e+2f,
        8.25831750e+5f, 2.66378633e+4f, 3.29637512e+2f, 8.80438438e+5f, 2.79782715e+4f, 3.40533386e+2f,
        9.37769688e+5f, 2.93625469e+4f, 3.51570709e+2f, 9.97913812e+5f, 3.07912520e+4f, 3.62747925e+2f,
        1.06096012e+6f, 3.22649473e+4f, 3.74063690e+2f, 1.12699938e+6f, 3.37841797e+4f, 3.85516907e+2f,
        1.19612300e+6f, 3.53495039e+4f, 3.97106781e+2f, 1.26842388e+6f, 3.69614570e+4f, 4.08832764e+2f,
        1.34399562e+6f, 3.86205898e+4f, 4.20694519e+2f, 1.42293325e+6f, 4.03274414e+4f, 4.32692139e+2f,
        1.50533275e+6f, 4.20825586e+4f, 4.44826019e+2f, 1.59129112e+6f, 4.38864844e+4f, 4.57096954e+2f,
        1.68090662e+6f, 4.57397734e+4f, 4.69506073e+2f, 1.77427850e+6f, 4.76429766e+4f, 4.82055023e+2f,
        1.87150712e+6f, 4.95966602e+4f, 4.94745880e+2f, 1.97269412e+6f, 5.16013984e+4f, 5.07581268e+2f,
        2.07794200e+6f, 5.36577773e+4f, 5.20564392e+2f, 2.18735475e+6f, 5.57663867e+4f, 5.33698975e+2f,
        2.30103750e+6f, 5.79278477e+4f, 5.46989624e+2f, 2.41909650e+6f, 6.01427969e+4f, 5.60441467e+2f,
        2.54163925e+6f, 6.24118828e+4f, 5.74060486e+2f, 2.66877500e+6f, 6.47357969e+4f, 5.87853760e+2f,
        2.80061400e+6f, 6.71152422e+4f, 6.01829041e+2f, 2.93726800e+6f, 6.95509766e+4f, 6.15995422e+2f,
        3.07885025e+6f, 7.20437734e+4f, 6.30362976e+2f, 3.22547575e+6f, 7.45944609e+4f, 6.44943359e+2f,
        3.37726150e+6f, 7.72039297e+4f, 6.59749451e+2f, 3.53432550e+6f, 7.98730938e+4f, 6.74795898e+2f,
        3.69678825e+6f, 8.26029531e+4f, 6.90099243e+2f, 3.86477225e+6f, 8.53945781e+4f, 7.05678101e+2f,
        4.03840200e+6f, 8.82491016e+4f, 7.21553162e+2f, 4.21780500e+6f, 9.11677656e+4f, 7.37748169e+2f,
        4.40311050e+6f, 9.41518984e+4f, 7.54289795e+2f, 4.59445050e+6f, 9.72029453e+4f, 7.71208130e+2f,
        4.79196100e+6f, 1.00322477e+5f, 7.88537292e+2f, 4.99578000e+6f, 1.03512219e+5f, 8.06316345e+2f,
        5.20605050e+6f, 1.06774055e+5f, 8.24590088e+2f, 5.42291850e+6f, 1.10110070e+5f, 8.43410278e+2f,
        5.64653400e+6f, 1.13522562e+5f, 8.62837280e+2f, 5.87705350e+6f, 1.17014109e+5f, 8.82942017e+2f,
        6.11463700e+6f, 1.20587570e+5f, 9.03808960e+2f, 6.35945200e+6f, 1.24246203e+5f, 9.25539246e+2f,
        6.61167200e+6f, 1.27993695e+5f, 9.48255249e+2f, 6.87147950e+6f, 1.31834281e+5f, 9.72106384e+2f,
        7.13906500e+6f, 1.35772875e+5f, 9.97277039e+2f, 7.41463000e+6f, 1.39815172e+5f, 1.02399915e+3f,
        7.69838800e+6f, 1.43967984e+5f, 1.05256982e+3f, 7.99056900e+6f, 1.48239438e+5f, 1.08338110e+3f,
        8.29141900e+6f, 1.52639531e+5f, 1.11696545e+3f, 8.60120700e+6f, 1.57180750e+5f, 1.15407227e+3f,
        8.92023100e+6f, 1.61879281e+5f, 1.19579395e+3f, 9.24882500e+6f, 1.66756688e+5f, 1.24379053e+3f,
        9.58737600e+6f, 1.71843000e+5f, 1.30070911e+3f, 9.93634000e+6f, 1.77182156e+5f, 1.37105212e+3f,
        1.02962850e+7f, 1.82843047e+5f, 1.46321106e+3f, 1.06679580e+7f, 1.88944641e+5f, 1.59559485e+3f,
    };

    private static ReadOnlySpan<float> CoeffHigh => new float[]
    {
        1.08220520e+7f, 4.78906758e+4f, 1.05626465e+2f, 1.09182560e+7f, 4.83157070e+4f, 1.07561577e+2f,
        1.10153230e+7f, 4.87538789e+4f, 1.11550323e+2f, 1.11132820e+7f, 4.92096250e+4f, 1.16380981e+2f,
        1.12121730e+7f, 4.96872227e+4f, 1.22578972e+2f, 1.13135560e+7f, 3.83243359e+4f, -1.86568145e+4f,
    };
}
