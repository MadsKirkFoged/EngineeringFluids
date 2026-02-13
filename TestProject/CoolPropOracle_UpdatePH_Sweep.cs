//using Microsoft.VisualStudio.TestTools.UnitTesting;
//using EngineeringUnits;
//using System;
//using System.Collections.Generic;
//using System.Globalization;
//using System.IO;
//using EngineeringFluids.Helmholtz.Solvers;
//using EngineeringFluids.Fluids;

//[TestClass]
//public class CoolPropOracle_UpdatePH_Sweep
//{
//    private static double RelErr(double expected, double actual)
//    {
//        double denom = Math.Max(Math.Abs(expected), 1.0);
//        return Math.Abs(actual - expected) / denom;
//    }

//    private static double Clamp(double x, double lo, double hi) => x < lo ? lo : (x > hi ? hi : x);

//    private static double LogUniform(Random rng, double a, double b)
//    {
//        double la = Math.Log(a);
//        double lb = Math.Log(b);
//        return Math.Exp(la + (lb - la) * rng.NextDouble());
//    }

//    private static void SkipIfOracleFailed(SharpFluids.Fluid f, string context)
//    {
//        if (f.FailState)
//            Assert.Inconclusive($"SharpFluids oracle failed (FailState=true): {context}");
//    }

//    private static double Percentile(double[] data, double p)
//    {
//        if (data.Length == 0)
//            return double.NaN;
//        Array.Sort(data);
//        double pos = (data.Length - 1) * p;
//        int i = (int)pos;
//        double frac = pos - i;
//        if (i + 1 < data.Length)
//            return data[i] * (1 - frac) + data[i + 1] * frac;
//        return data[i];
//    }

//    [TestMethod]
//    [TestCategory("LongRunning")]
//    public void Sweep_UpdatePH_1000_RoundTrip_EdgeAware()
//    {
//        // =========================
//        // Config
//        // =========================
//        const int N = 1000; // set to 100_000 for big run

//        // "Normal" single-phase tolerances (gating)
//        const double relTolProps_1ph = 1e-5;
//        const double absTolP_Pa_1ph = 5.0;
//        const double absTolT_K_1ph = 5e-3;
//        const double relTolRoundTripH_1ph = 1e-8;

//        // Two-phase tolerances (gating) – do NOT use PT roundtrip in 2ph (ill-posed)
//        const double relTolH_2ph = 1e-6;
//        const double absTolT_K_2ph = 2e-2;   // 0.02 K
//        const double absTolP_Pa_2ph = 20.0;  // 20 Pa

//        // (A) Stop gating on q match; just require valid range:
//        const bool gateOnQualityMatch_2ph = false;
//        const double absTolQ_2ph = 1e-3; // only used if gateOnQualityMatch_2ph == true

//        // (B) Stop gating on rho in two-phase:
//        const bool gateOnRho_2ph = false;
//        const double relTolRho_2ph = 5e-3; // only used if gateOnRho_2ph == true

//        // Narrow skip bands OR convert to bucket-based tolerances
//        const double nearSatBand_K = 0.01;
//        const double nearCriticalBand_K = 0.2;
//        const double nearCriticalBand_P = 0.995;


//        const bool writeCsv = true;
//        const int keepWorstN = 200;

//        // Our fluid
//        var a = new EngineeringFluids.Fluids.Ammonia();

//        // Sampling bounds
//        double PcPa = new Ammonia().Critical.Pressure.Pascal;
//        double TcK = new Ammonia().Critical.Temperature.Kelvin;

//        const double PminPa = 5_000.0;
//        double PsubMax = 0.999 * PcPa;
//        double PsuperMin = 1.05 * PcPa;
//        double PsuperMax = 3.0 * PcPa;

//        const double TminK = 200.0;
//        const double TmaxK = 1200.0;

//        // ΔT distribution around Tsat(P)
//        const double dT_near_min = 1e-3;
//        const double dT_near_max = 1.0;
//        const double dT_far_min = 1.0;
//        const double dT_far_max = 200.0;

//        // Two-phase quality range
//        const double qMin = 0.05;
//        const double qMax = 0.95;

//        // =========================
//        // Objects
//        // =========================
//        var rng = new Random(1234567);

//        // Oracle fluids (CoolProp via SharpFluids)
//        var refFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        var satFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        var L = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        var V = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);



//        // Stats
//        var eRho = new List<double>(N);
//        var eH = new List<double>(N);
//        var eS = new List<double>(N);
//        var eU = new List<double>(N);
//        var eT = new List<double>(N);
//        var eP = new List<double>(N);
//        var eHRt = new List<double>(N); // single-phase only

//        int pass = 0, fail = 0, skip = 0;
//        int skipNearSat = 0, skipNearCritical = 0;

//        // Worst points
//        var worst = new List<(double eRho, string region, double P, double Tref, double Tcalc,
//                              double rhoRef, double rho, double hRef, double h,
//                              double qRef, double q, double eHRt, string reason)>(keepWorstN);

//        // CSV
//        List<string> lines = null;
//        if (writeCsv)
//        {
//            lines = new List<string>(N + 1);
//            lines.Add("region,P_Pa,Tref_K,Tcalc_K,rho_ref,rho,err_rho,h_ref,h,err_h,s_ref,s,err_s,u_ref,u,err_u,p_ref,p,err_pabs,q_ref,q,err_h_roundtrip,status,reason");
//        }

//        // =========================
//        // Sampling helpers
//        // =========================
//        double SampleSubcriticalPressure()
//        {
//            double u = rng.NextDouble();
//            if (u < 0.60)
//            {
//                return Math.Exp(Math.Log(PminPa) + rng.NextDouble() * (Math.Log(0.60 * PcPa) - Math.Log(PminPa)));
//            }
//            else
//            {
//                double v = rng.NextDouble();
//                v = 1.0 - Math.Pow(1.0 - v, 3.0);
//                return (0.60 * PcPa) + v * (PsubMax - 0.60 * PcPa);
//            }
//        }

//        double SampleSupercriticalPressure()
//        {
//            return Math.Exp(Math.Log(PsuperMin) + rng.NextDouble() * (Math.Log(PsuperMax) - Math.Log(PsuperMin)));
//        }

//        double SampleSupercriticalTemperature()
//        {
//            double lo = TcK + 1.0;
//            double hi = Math.Min(TmaxK, TcK + 500.0);
//            return lo + rng.NextDouble() * (hi - lo);
//        }

//        // =========================
//        // Main loop
//        // =========================
//        int idx = 0;
//        while (idx < N)
//        {
//            double r = rng.NextDouble();

//            string region;
//            double PPa;
//            double TrefK = double.NaN;
//            double qRef = -1.0;

//            Pressure P;
//            Temperature Tref;

//            // Oracle reference values
//            double rho_ref, h_ref, s_ref, u_ref, p_ref, T_ref_out;

//            // For edge checks
//            double TsatOracleK = double.NaN;
//            bool haveTsatOracle = false;

//            try
//            {
//                if (r < 0.70)
//                {
//                    // Subcritical single-phase around Tsat(P)
//                    region = (rng.NextDouble() < 0.5) ? "vapor" : "liquid";
//                    PPa = SampleSubcriticalPressure();
//                    P = Pressure.FromPascal(PPa);

//                    // Oracle Tsat(P)
//                    satFluid.Pressure = P;
//                    var TsatObj = satFluid.Tsat;
//                    SkipIfOracleFailed(satFluid, $"Tsat at P={PPa} Pa");
//                    if (TsatObj == null)
//                    { skip++; continue; }

//                    TsatOracleK = TsatObj.Kelvin;
//                    haveTsatOracle = true;

//                    double dT = (rng.NextDouble() < 0.5)
//                        ? LogUniform(rng, dT_near_min, dT_near_max)
//                        : LogUniform(rng, dT_far_min, dT_far_max);

//                    double TK = region == "vapor" ? (TsatOracleK + dT) : (TsatOracleK - dT);
//                    TK = Clamp(TK, TminK, TmaxK);

//                    TrefK = TK;
//                    Tref = Temperature.FromKelvin(TK);

//                    refFluid.UpdatePT(P, Tref);
//                    SkipIfOracleFailed(refFluid, $"UpdatePT oracle at P={PPa} Pa, T={TK} K");

//                    // If oracle ended up two-phase (tiny ΔT), skip
//                    if (refFluid.Quality >= 0.0 && refFluid.Quality <= 1.0)
//                    {
//                        skip++;
//                        continue;
//                    }

//                    rho_ref = refFluid.Density!.KilogramPerCubicMeter;
//                    h_ref = refFluid.Enthalpy!.JoulePerKilogram;
//                    s_ref = refFluid.Entropy!.JoulePerKilogramKelvin;
//                    u_ref = refFluid.InternalEnergy!.JoulePerKilogram;
//                    p_ref = refFluid.Pressure!.Pascal;
//                    T_ref_out = refFluid.Temperature!.Kelvin;
//                    qRef = -1.0;
//                }
//                else if (r < 0.85)
//                {
//                    // Two-phase
//                    region = "twophase";
//                    PPa = SampleSubcriticalPressure();
//                    PPa = Math.Min(PPa, 0.95 * PcPa);
//                    P = Pressure.FromPascal(PPa);

//                    qRef = qMin + rng.NextDouble() * (qMax - qMin);

//                    // Endpoints
//                    L.UpdatePX(P, 0.0);
//                    SkipIfOracleFailed(L, $"UpdatePX(P,0) oracle at P={PPa} Pa");
//                    V.UpdatePX(P, 1.0);
//                    SkipIfOracleFailed(V, $"UpdatePX(P,1) oracle at P={PPa} Pa");

//                    if (L.Enthalpy == null || V.Enthalpy == null)
//                    { skip++; continue; }

//                    double hL_ref = L.Enthalpy.JoulePerKilogram;
//                    double hV_ref = V.Enthalpy.JoulePerKilogram;

//                    if (Math.Abs(hV_ref - hL_ref) < 1e-2)
//                    { skip++; continue; } // ill-conditioned

//                    double hMix = (1.0 - qRef) * hL_ref + qRef * hV_ref;

//                    refFluid.UpdatePX(P, qRef);
//                    SkipIfOracleFailed(refFluid, $"UpdatePX oracle at P={PPa} Pa, q={qRef}");

//                    rho_ref = refFluid.Density!.KilogramPerCubicMeter;
//                    h_ref = hMix;
//                    s_ref = refFluid.Entropy!.JoulePerKilogramKelvin;
//                    u_ref = refFluid.InternalEnergy!.JoulePerKilogram;
//                    p_ref = refFluid.Pressure!.Pascal;
//                    T_ref_out = refFluid.Temperature!.Kelvin;

//                    // Tsat is basically T in two-phase oracle state
//                    TsatOracleK = T_ref_out;
//                    haveTsatOracle = true;

//                    TrefK = T_ref_out;
//                    Tref = Temperature.FromKelvin(TrefK);
//                }
//                else
//                {
//                    // Supercritical
//                    region = "supercritical";
//                    PPa = SampleSupercriticalPressure();
//                    P = Pressure.FromPascal(PPa);

//                    double TK = SampleSupercriticalTemperature();
//                    TK = Clamp(TK, TminK, TmaxK);

//                    TrefK = TK;
//                    Tref = Temperature.FromKelvin(TK);

//                    refFluid.UpdatePT(P, Tref);
//                    SkipIfOracleFailed(refFluid, $"UpdatePT oracle at P={PPa} Pa, T={TK} K");

//                    rho_ref = refFluid.Density!.KilogramPerCubicMeter;
//                    h_ref = refFluid.Enthalpy!.JoulePerKilogram;
//                    s_ref = refFluid.Entropy!.JoulePerKilogramKelvin;
//                    u_ref = refFluid.InternalEnergy!.JoulePerKilogram;
//                    p_ref = refFluid.Pressure!.Pascal;
//                    T_ref_out = refFluid.Temperature!.Kelvin;
//                    qRef = -1.0;
//                }
//            }
//            catch
//            {
//                skip++;
//                continue;
//            }

//            // Our solve
//            try
//            {
//                a.UpdatePH(P, Enthalpy.FromJoulePerKilogram(h_ref));

//                double rho = a.Density!.KilogramPerCubicMeter;
//                double h = a.Enthalpy.SI;
//                double s = a.Entropy.SI;
//                double u = a.InternalEnergy.SI;
//                double pCalc = a.Pressure.Pascal;
//                double Tcalc = a.Temperature!.Kelvin;
//                double q = a.Quality;

//                double erho = RelErr(rho_ref, rho);
//                double eh = RelErr(h_ref, h);
//                double es = RelErr(s_ref, s);
//                double eu = RelErr(u_ref, u);
//                double eTloc = Math.Abs(Tcalc - T_ref_out);
//                double ep = Math.Abs(pCalc - p_ref);

//                bool isTwoPhaseCase = (region == "twophase");
//                bool nearCritical = (P.Pascal > nearCriticalBand_P * PcPa) || (Math.Abs(T_ref_out - TcK) < nearCriticalBand_K);

//                bool nearSat = false;
//                if (!isTwoPhaseCase && haveTsatOracle)
//                {
//                    nearSat = Math.Abs(T_ref_out - TsatOracleK) < nearSatBand_K;
//                }

//                // If oracle says single-phase but we're too near critical/sat, we *log* but don't fail the sweeper.
//                // This avoids false negatives in knife-edge regions where flash routines are numerically delicate. [1](https://github.com/CoolProp/CoolProp/discussions/2403)[2](https://colab.research.google.com/github/JMQuinlan/Thermo/blob/main/Thermo_Calc_Example.ipynb)
//                if (!isTwoPhaseCase && (nearCritical || nearSat))
//                {
//                    skip++;
//                    if (nearCritical)
//                        skipNearCritical++;
//                    if (nearSat)
//                        skipNearSat++;

//                    if (writeCsv)
//                    {
//                        lines.Add(string.Join(",",
//                            region,
//                            P.Pascal.ToString("G17", CultureInfo.InvariantCulture),
//                            T_ref_out.ToString("G17", CultureInfo.InvariantCulture),
//                            Tcalc.ToString("G17", CultureInfo.InvariantCulture),
//                            rho_ref.ToString("G17", CultureInfo.InvariantCulture),
//                            rho.ToString("G17", CultureInfo.InvariantCulture),
//                            erho.ToString("G17", CultureInfo.InvariantCulture),
//                            h_ref.ToString("G17", CultureInfo.InvariantCulture),
//                            h.ToString("G17", CultureInfo.InvariantCulture),
//                            eh.ToString("G17", CultureInfo.InvariantCulture),
//                            s_ref.ToString("G17", CultureInfo.InvariantCulture),
//                            s.ToString("G17", CultureInfo.InvariantCulture),
//                            es.ToString("G17", CultureInfo.InvariantCulture),
//                            u_ref.ToString("G17", CultureInfo.InvariantCulture),
//                            u.ToString("G17", CultureInfo.InvariantCulture),
//                            eu.ToString("G17", CultureInfo.InvariantCulture),
//                            p_ref.ToString("G17", CultureInfo.InvariantCulture),
//                            pCalc.ToString("G17", CultureInfo.InvariantCulture),
//                            ep.ToString("G17", CultureInfo.InvariantCulture),
//                            qRef.ToString("G17", CultureInfo.InvariantCulture),
//                            q.ToString("G17", CultureInfo.InvariantCulture),
//                            "",
//                            "SKIP",
//                            $"edge_skip nearSat={nearSat} nearCritical={nearCritical}"
//                        ));
//                    }

//                    idx++;
//                    continue;
//                }

//                bool ok;
//                double e_rt = double.NaN;

//                if (!isTwoPhaseCase)
//                {
//                    // Single-phase strict grading + PT roundtrip
//                    a.UpdatePT(P, Temperature.FromKelvin(Tcalc));
//                    double h_rt = a.Enthalpy.SI;
//                    e_rt = RelErr(h_ref, h_rt);

//                    ok =
//                        erho <= relTolProps_1ph &&
//                        eh <= relTolProps_1ph &&
//                        es <= relTolProps_1ph &&
//                        eu <= relTolProps_1ph &&
//                        ep <= absTolP_Pa_1ph &&
//                        eTloc <= absTolT_K_1ph &&
//                        e_rt <= relTolRoundTripH_1ph &&
//                        (q == -1.0);

//                    eHRt.Add(e_rt);
//                }
//                else
//                {
//                    // Two-phase grading: check q, Tsat, p, h
//                    ok =
//                        (q >= 0.0 && q <= 1.0) &&
//                        eh <= relTolH_2ph &&
//                        eTloc <= absTolT_K_2ph &&
//                        ep <= absTolP_Pa_2ph;

//                    if (gateOnQualityMatch_2ph)
//                        ok = ok && (Math.Abs(q - qRef) < absTolQ_2ph);

//                    if (gateOnRho_2ph)
//                        ok = ok && (erho <= relTolRho_2ph);

//                    // We DO NOT gate on rho/s/u tightly in twophase – endpoint dependent. [1](https://github.com/CoolProp/CoolProp/discussions/2403)[2](https://colab.research.google.com/github/JMQuinlan/Thermo/blob/main/Thermo_Calc_Example.ipynb)
//                }

//                if (ok)
//                    pass++;
//                else
//                    fail++;

//                // stats
//                eRho.Add(erho);
//                eH.Add(eh);
//                eS.Add(es);
//                eU.Add(eu);
//                eT.Add(eTloc);
//                eP.Add(ep);

//                // worst
//                string reason = ok ? "PASS" : (isTwoPhaseCase ? "FAIL_2PH" : "FAIL_1PH");
//                if (worst.Count < keepWorstN)
//                {
//                    worst.Add((erho, region, P.Pascal, T_ref_out, Tcalc, rho_ref, rho, h_ref, h, qRef, q, e_rt, reason));
//                    worst.Sort((a1, a2) => a2.eRho.CompareTo(a1.eRho));
//                }
//                else if (erho > worst[worst.Count - 1].eRho)
//                {
//                    worst[worst.Count - 1] = (erho, region, P.Pascal, T_ref_out, Tcalc, rho_ref, rho, h_ref, h, qRef, q, e_rt, reason);
//                    worst.Sort((a1, a2) => a2.eRho.CompareTo(a1.eRho));
//                }

//                if (writeCsv)
//                {
//                    lines.Add(string.Join(",",
//                        region,
//                        P.Pascal.ToString("G17", CultureInfo.InvariantCulture),
//                        T_ref_out.ToString("G17", CultureInfo.InvariantCulture),
//                        Tcalc.ToString("G17", CultureInfo.InvariantCulture),
//                        rho_ref.ToString("G17", CultureInfo.InvariantCulture),
//                        rho.ToString("G17", CultureInfo.InvariantCulture),
//                        erho.ToString("G17", CultureInfo.InvariantCulture),
//                        h_ref.ToString("G17", CultureInfo.InvariantCulture),
//                        h.ToString("G17", CultureInfo.InvariantCulture),
//                        eh.ToString("G17", CultureInfo.InvariantCulture),
//                        s_ref.ToString("G17", CultureInfo.InvariantCulture),
//                        s.ToString("G17", CultureInfo.InvariantCulture),
//                        es.ToString("G17", CultureInfo.InvariantCulture),
//                        u_ref.ToString("G17", CultureInfo.InvariantCulture),
//                        u.ToString("G17", CultureInfo.InvariantCulture),
//                        eu.ToString("G17", CultureInfo.InvariantCulture),
//                        p_ref.ToString("G17", CultureInfo.InvariantCulture),
//                        pCalc.ToString("G17", CultureInfo.InvariantCulture),
//                        ep.ToString("G17", CultureInfo.InvariantCulture),
//                        qRef.ToString("G17", CultureInfo.InvariantCulture),
//                        q.ToString("G17", CultureInfo.InvariantCulture),
//                        (double.IsNaN(e_rt) ? "" : e_rt.ToString("G17", CultureInfo.InvariantCulture)),
//                        ok ? "PASS" : "FAIL",
//                        reason
//                    ));
//                }
//            }
//            catch (Exception ex)
//            {
//                fail++;
//                if (writeCsv)
//                {
//                    lines.Add($"{region},{P.Pascal.ToString(CultureInfo.InvariantCulture)},{TrefK.ToString(CultureInfo.InvariantCulture)},," +
//                              $"{rho_ref.ToString(CultureInfo.InvariantCulture)},,,{h_ref.ToString(CultureInfo.InvariantCulture)},,,{s_ref.ToString(CultureInfo.InvariantCulture)},,,{u_ref.ToString(CultureInfo.InvariantCulture)},," +
//                              $"{p_ref.ToString(CultureInfo.InvariantCulture)},,,{qRef.ToString(CultureInfo.InvariantCulture)},,," +
//                              $"FAIL,\"THROW: {ex.GetType().Name}: {ex.Message.Replace('\"', '\'')}\"");
//                }
//            }

//            idx++;
//        }

//        // =========================
//        // Summary
//        // =========================
//        string Summary(List<double> v, string name)
//        {
//            var arr = v.ToArray();
//            return $"{name} p50/p95/p99/p999 = {Percentile(arr, 0.50):g3} / {Percentile(arr, 0.95):g3} / {Percentile(arr, 0.99):g3} / {Percentile(arr, 0.999):g3}";
//        }

//        string summary =
//            $"PH Sweep done: points={N}, pass={pass}, fail={fail}, skip={skip} (skipNearSat={skipNearSat}, skipNearCritical={skipNearCritical})\n" +
//            Summary(eRho, "err_rho") + "\n" +
//            Summary(eH, "err_h") + "\n" +
//            Summary(eS, "err_s") + "\n" +
//            Summary(eU, "err_u") + "\n" +
//            Summary(eT, "abs_err_T(K)") + "\n" +
//            Summary(eP, "abs_err_p(Pa)") + "\n";

//        if (eHRt.Count > 0)
//            summary += Summary(eHRt, "roundtrip_err_h (single-phase only)") + "\n";

//        summary += "\nWorst points (by rho rel error):\n";
//        foreach (var w in worst)
//        {
//            summary += $"  {w.region} P={w.P:g0} Pa, Tref={w.Tref:g6} K, T={w.Tcalc:g6} K, " +
//                       $"rho_ref={w.rhoRef:g6}, rho={w.rho:g6}, err_rho={w.eRho:g3}, " +
//                       $"h_ref={w.hRef:g6}, h={w.h:g6}, q_ref={w.qRef:g3}, q={w.q:g3}, rt_h_err={w.eHRt:g3}, reason={w.reason}\n";
//        }

//        if (writeCsv)
//        {
//            string csvPath = Path.Combine(Directory.GetCurrentDirectory(), $"PH_sweep_{N}_{DateTime.UtcNow:yyyyMMdd_HHmmss}.csv");
//            File.WriteAllLines(csvPath, lines);
//            summary += $"\nCSV written: {csvPath}\n";
//        }

//        if (fail > 0)
//            Assert.Fail(summary);
//    }


//    [TestMethod]
//    [TestCategory("LongRunning")]
//    public void Sweep_UpdatePX_5000_Oracle_Diagnostic_v2()
//    {
//        const int N = 5000;

//        // Diagnostic mode: always completes and prints summary
//        const bool diagnosticOnly = true;

//        // Two-phase tolerances vs CoolProp oracle
//        const double absTolP_Pa = 2.0;       // should be tight (PQ flash enforces P)
//        const double absTolT_K = 5e-2;      // 0.05 K
//        const double relTolH = 1e-3;      // 0.01%
//        const double relTolRho = 5e-3;      // 0.5%
//        const double absTolQ = 1e-12;     // q is set exactly on our side

//        double PcPa = new Ammonia().Critical.Pressure.Pascal;

//        var rng = new Random(1234567);
//        var refFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);

//        int pass = 0, fail = 0, skip = 0;
//        int failP = 0, failT = 0, failH = 0, failRho = 0, failQ = 0;

//        var examples = new List<string>(20);

//        var lines = new List<string>(N + 1);
//        lines.Add("P_Pa,q,Tsat_ref_K,Tcalc_K,absErrT_K,h_ref,h,relErrH,rho_ref,rho,relErrRho,p_ref,p,absErrP,status,reason");

//        for (int i = 0; i < N; i++)
//        {
//            // sample pressure below critical
//            double PPa = Math.Exp(Math.Log(5000.0) + rng.NextDouble() * (Math.Log(0.95 * PcPa) - Math.Log(5000.0)));
//            var P = Pressure.FromPascal(PPa);

//            // random quality away from endpoints
//            double q = 0.05 + rng.NextDouble() * 0.90;

//            // Oracle state at (P,q) – use its Temperature as Tsat reference
//            refFluid.UpdatePX(P, q);
//            if (refFluid.FailState || refFluid.Enthalpy == null || refFluid.Density == null || refFluid.Temperature == null || refFluid.Pressure == null)
//            {
//                skip++;
//                continue;
//            }

//            double Tsat_ref = refFluid.Temperature.Kelvin;
//            double h_ref = refFluid.Enthalpy.JoulePerKilogram;
//            double rho_ref = refFluid.Density.KilogramPerCubicMeter;
//            double p_ref = refFluid.Pressure.Pascal;

//            // Our UpdatePX
//            var a = new EngineeringFluids.Fluids.Ammonia();
//            try
//            {
//                a.UpdatePX(P, q);
//            }
//            catch
//            {
//                fail++;
//                lines.Add($"{PPa.ToString(CultureInfo.InvariantCulture)},{q.ToString(CultureInfo.InvariantCulture)},{Tsat_ref.ToString(CultureInfo.InvariantCulture)},,,,,,,,,,,FAIL,THREW");
//                continue;
//            }

//            double T = a.Temperature!.Kelvin;
//            double p = a.Pressure.Pascal;
//            double h = a.Enthalpy.SI;
//            double rho = a.Density!.KilogramPerCubicMeter;
//            double qCalc = a.Quality;

//            double absErrT = Math.Abs(T - Tsat_ref);
//            double absErrP = Math.Abs(p - PPa);
//            double relErrH = RelErr(h_ref, h);
//            double relErrRho = RelErr(rho_ref, rho);
//            double absErrQ = Math.Abs(qCalc - q);

//            bool ok = true;
//            string reason = "";

//            if (absErrP > absTolP_Pa)
//            { ok = false; reason += "P "; failP++; }
//            if (absErrT > absTolT_K)
//            { ok = false; reason += "T "; failT++; }
//            if (relErrH > relTolH)
//            { ok = false; reason += "H "; failH++; }
//            if (relErrRho > relTolRho)
//            { ok = false; reason += "RHO "; failRho++; }
//            if (absErrQ > absTolQ)
//            { ok = false; reason += "Q "; failQ++; }

//            if (ok)
//                pass++;
//            else
//                fail++;

//            if (!ok && examples.Count < 10)
//            {
//                examples.Add($"FAIL P={PPa:g0} Pa q={q:g3} | dT={absErrT:g3}K relH={relErrH:g3} relRho={relErrRho:g3} dP={absErrP:g3}Pa qErr={absErrQ:g3} [{reason.Trim()}]");
//            }

//            lines.Add(string.Join(",",
//                PPa.ToString("G17", CultureInfo.InvariantCulture),
//                q.ToString("G17", CultureInfo.InvariantCulture),
//                Tsat_ref.ToString("G17", CultureInfo.InvariantCulture),
//                T.ToString("G17", CultureInfo.InvariantCulture),
//                absErrT.ToString("G17", CultureInfo.InvariantCulture),
//                h_ref.ToString("G17", CultureInfo.InvariantCulture),
//                h.ToString("G17", CultureInfo.InvariantCulture),
//                relErrH.ToString("G17", CultureInfo.InvariantCulture),
//                rho_ref.ToString("G17", CultureInfo.InvariantCulture),
//                rho.ToString("G17", CultureInfo.InvariantCulture),
//                relErrRho.ToString("G17", CultureInfo.InvariantCulture),
//                p_ref.ToString("G17", CultureInfo.InvariantCulture),
//                p.ToString("G17", CultureInfo.InvariantCulture),
//                absErrP.ToString("G17", CultureInfo.InvariantCulture),
//                ok ? "PASS" : "FAIL",
//                $"\"{reason.Trim()}\""
//            ));
//        }

//        string csvPath = Path.Combine(Directory.GetCurrentDirectory(), $"PX_sweep_{N}_{DateTime.UtcNow:yyyyMMdd_HHmmss}.csv");
//        File.WriteAllLines(csvPath, lines);

//        string summary =
//            $"UpdatePX sweep: N={N}, pass={pass}, fail={fail}, skip={skip}\n" +
//            $"Fail breakdown: P={failP}, T={failT}, H={failH}, RHO={failRho}, Q={failQ}\n" +
//            $"CSV: {csvPath}\n" +
//            (examples.Count > 0 ? ("Examples:\n  " + string.Join("\n  ", examples)) : "");

//        if (!diagnosticOnly && fail > 0)
//            Assert.Fail(summary);

//        Assert.Inconclusive(summary);
//    }



//    //private static double RelErr(double expected, double actual)
//    //{
//    //    double denom = Math.Max(Math.Abs(expected), 1.0);
//    //    return Math.Abs(actual - expected) / denom;
//    //}

//    //private static void SkipIfOracleFailed(SharpFluids.Fluid f, string context)
//    //{
//    //    if (f.FailState)
//    //        Assert.Inconclusive($"SharpFluids oracle failed (FailState=true): {context}");
//    //}

//    // -----------------------------
//    // Single-phase test points
//    // -----------------------------
//    public static IEnumerable<object[]> SinglePhasePoints()
//    {
//        // (Pbar, mode, deltaT_K) where mode is "vapor"/"liquid"
//        double[] pBars = { 1, 5, 10, 20, 40, 60, 80, 100 };
//        double[] dTs = { 2, 10, 30, 60 };

//        foreach (var p in pBars)
//        {
//            foreach (var dT in dTs)
//            {
//                yield return new object[] { p, "vapor", dT };
//                yield return new object[] { p, "liquid", dT };
//            }
//        }

//        // Supercritical absolute points (Pbar, mode="abs", T_K)
//        yield return new object[] { 150.0, "abs", 450.0 };
//        yield return new object[] { 250.0, "abs", 600.0 };
//        yield return new object[] { 300.0, "abs", 750.0 };
//    }

//    [DataTestMethod]
//    [DynamicData(nameof(SinglePhasePoints), DynamicDataSourceType.Method)]
//    public void UpdatePS_SinglePhase_Matches_SharpFluids_AndQuality(double pBar, string mode, double value)
//    {
//        // Tolerances similar to PH tests; adjust if you want tighter later
//        const double relTolProps = 1e-5;
//        const double absTolP_Pa = 5.0;
//        const double absTolT_K = 5e-3;

//        Pressure P = Pressure.FromBar(pBar);
//        Temperature T;

//        // Oracle Tsat(P) for constructing vapor/liquid points
//        if (mode == "abs")
//        {
//            T = Temperature.FromKelvin(value);
//        }
//        else
//        {
//            var satFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//            satFluid.Pressure = P;
//            var TsatObj = satFluid.Tsat;
//            SkipIfOracleFailed(satFluid, $"Tsat at P={P.Pascal} Pa");
//            if (TsatObj == null)
//                Assert.Inconclusive($"Tsat null at P={P.Pascal} Pa");

//            double TsatK = TsatObj!.Kelvin;
//            double dT = value;

//            double TK = mode == "vapor" ? (TsatK + dT) : (TsatK - dT);
//            if (TK < 200.0)
//                Assert.Inconclusive($"Skipping too-low T={TK} K");
//            T = Temperature.FromKelvin(TK);
//        }

//        // Oracle reference at (P,T)
//        var refFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        refFluid.UpdatePT(P, T);
//        SkipIfOracleFailed(refFluid, $"UpdatePT oracle at P={P.Pascal} Pa, T={T.Kelvin} K");

//        // Skip if oracle landed in two-phase
//        if (refFluid.Quality >= 0.0 && refFluid.Quality <= 1.0)
//            Assert.Inconclusive($"Oracle two-phase at P={P.Pascal} Pa, T={T.Kelvin} K (Q={refFluid.Quality})");

//        double s_ref = refFluid.Entropy!.JoulePerKilogramKelvin;
//        double rho_ref = refFluid.Density!.KilogramPerCubicMeter;
//        double h_ref = refFluid.Enthalpy!.JoulePerKilogram;
//        double u_ref = refFluid.InternalEnergy!.JoulePerKilogram;
//        double p_ref = refFluid.Pressure!.Pascal;
//        double T_ref = refFluid.Temperature!.Kelvin;

//        // Your flash
//        var a = new EngineeringFluids.Fluids.Ammonia();
//        a.UpdatePS(P, SpecificEntropy.FromJoulePerKilogramKelvin(s_ref));

//        // Single-phase contract
//        Assert.IsTrue(a.Quality == -1.0, $"Expected Quality=-1 in single-phase; got {a.Quality}");

//        // Compare
//        Assert.IsTrue(RelErr(rho_ref, a.Density!.KilogramPerCubicMeter) <= relTolProps, "rho mismatch");
//        Assert.IsTrue(RelErr(h_ref, a.Enthalpy.SI) <= relTolProps, "h mismatch");
//        Assert.IsTrue(RelErr(u_ref, a.InternalEnergy.SI) <= relTolProps, "u mismatch");
//        Assert.IsTrue(RelErr(s_ref, a.Entropy.SI) <= relTolProps, "s mismatch");

//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - p_ref) < absTolP_Pa, "p mismatch vs oracle");
//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - P.Pascal) < absTolP_Pa, "p mismatch vs target");
//        Assert.IsTrue(Math.Abs(a.Temperature!.Kelvin - T_ref) < absTolT_K, "T mismatch");
//    }

//    // -----------------------------
//    // Two-phase test points
//    // -----------------------------
//    public static IEnumerable<object[]> TwoPhasePoints_Safe()
//    {
//        // Keep below critical margin; avoid endpoints
//        double[] pBars = { 1, 2, 5, 10, 20, 40, 60, 80 };
//        double[] qs = { 0.2, 0.5, 0.8 };

//        foreach (var p in pBars)
//            foreach (var q in qs)
//                yield return new object[] { p, q };
//    }

//    [DataTestMethod]
//    [DynamicData(nameof(TwoPhasePoints_Safe), DynamicDataSourceType.Method)]
//    public void UpdatePS_TwoPhase_Matches_SharpFluids_AndQuality(double pBar, double q)
//    {
//        // Two-phase tolerances: entropy and quality should match well; T and P are saturation constraints
//        const double absTolP_Pa = 20.0;
//        const double absTolT_K = 2e-2;
//        const double relTolS = 1e-6;
//        const double absTolQ = 1e-3;

//        Pressure P = Pressure.FromBar(pBar);

//        // Oracle endpoints
//        var L = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        var V = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);

//        L.UpdatePX(P, 0.0);
//        SkipIfOracleFailed(L, $"UpdatePX(P,0) oracle at P={P.Pascal} Pa");

//        V.UpdatePX(P, 1.0);
//        SkipIfOracleFailed(V, $"UpdatePX(P,1) oracle at P={P.Pascal} Pa");

//        if (L.Entropy == null || V.Entropy == null || L.Enthalpy == null || V.Enthalpy == null)
//            Assert.Inconclusive("Oracle returned null endpoints.");

//        double sL = L.Entropy.JoulePerKilogramKelvin;
//        double sV = V.Entropy.JoulePerKilogramKelvin;

//        // If ds collapses near critical, skip (ill-conditioned quality)
//        if (Math.Abs(sV - sL) < 1e-6)
//            Assert.Inconclusive("Skipping near-critical: sV≈sL");

//        double sMix = (1.0 - q) * sL + q * sV;

//        // Oracle two-phase reference
//        var refFluid = new SharpFluids.Fluid(SharpFluids.FluidList.Ammonia);
//        refFluid.UpdatePX(P, q);
//        SkipIfOracleFailed(refFluid, $"UpdatePX(P,q) oracle at P={P.Pascal} Pa, q={q}");

//        double Tsat_ref = refFluid.Temperature!.Kelvin;
//        double p_ref = refFluid.Pressure!.Pascal;

//        // Your flash
//        var a = new EngineeringFluids.Fluids.Ammonia();
//        a.UpdatePS(P, SpecificEntropy.FromJoulePerKilogramKelvin(sMix));

//        // Two-phase contract
//        Assert.IsTrue(a.Quality >= 0.0 && a.Quality <= 1.0, $"Expected two-phase Quality in [0,1], got {a.Quality}");
//        Assert.IsTrue(Math.Abs(a.Quality - q) < absTolQ, $"Quality mismatch: expected={q}, actual={a.Quality}");

//        // Saturation constraints
//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - P.Pascal) < absTolP_Pa, "p mismatch vs target");
//        Assert.IsTrue(Math.Abs(a.Pressure.Pascal - p_ref) < absTolP_Pa, "p mismatch vs oracle");
//        Assert.IsTrue(Math.Abs(a.Temperature!.Kelvin - Tsat_ref) < absTolT_K, "T mismatch vs Tsat");

//        // Entropy invariant
//        Assert.IsTrue(RelErr(sMix, a.Entropy.SI) < relTolS, $"s mismatch: expected={sMix}, actual={a.Entropy.SI}");
//    }

//}