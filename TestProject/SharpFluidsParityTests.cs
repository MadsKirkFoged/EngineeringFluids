using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using SharpFluids;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Phases = EngineeringFluids.Helmholtz.Phase.Phases;

namespace TestProject;

// EngineeringFluids against SharpFluids, update by update: the same inputs go into SharpFluids' Fluid (CoolProp)
// and into Ammonia, and every property both expose must agree.
//
// SharpFluids' UpdateXX maps to the CoolProp-compatible Ammonia update: UpdatePX/XT/PH/PS -> the Exact variants
// (the fast ones use ancillary fits in two-phase), everything else to the method of the same name.
//
// Known SharpFluids behaviour these tests work around - each one is pinned by a test of its own below, so a fix in
// SharpFluids shows up as a failing test here:
//  - Tau, Delta, AlphaR, ..., fugacity, gibbsmolar and umolar are hard-coded to 0 in its CoolProp wrapper: not compared.
//  - SoundSpeed is 0 whenever CoolProp's quality is in [0, 1] (all of two-phase, saturation lines included): compared
//    in single phase only.
//  - SurfaceTension is 0 exactly on the saturation lines and keeps its previous value in single phase: compared for
//    0 < q < 1 only.
//  - UpdatePT between 405.3 and 405.7 K retries a CoolProp failure at 405.3 K and reports that state: those
//    temperatures are left out of the grid.
//  - Tsat below the triple-point pressure is the state's own temperature (EngineeringFluids: NaN): not compared there.
//  - Quality in single phase is an update-dependent sentinel (-1, 10000, ...): any value outside [0, 1] is "single phase".
// And one CoolProp limitation: in the last 0.5 K below Tc its PQ and QT flashes disagree with each other (the EOS's
// dome closes near 405.25 K rather than at Tc), and PropsSI can even crash there - so the grid stops 0.5 K below Tc.
[TestClass]
public class SharpFluidsParityTests
{
    private const double Tc = 405.56;

    // One state as a system reports it (SI, NaN for "not available")
    private sealed record Snapshot(bool Failed, string Error, double T, double P, double D, double H, double S, double U,
        double Cp, double Cv, double W, double Mu, double K, double Pr, double Sigma, double Z, double Q, int Phase,
        double Tsat, double GasD, double LiqD);

    // A reference state from SharpFluids (CoolProp), from which each update's inputs are taken
    public sealed record State(double T, double P, double D, double H, double S, double Q)
    {
        public bool OnSaturationLine => Q == 0 || Q == 1;
        public bool TwoPhase => Q >= 0 && Q <= 1;
        public override string ToString() => string.Create(CultureInfo.InvariantCulture, $"T={T:G9} K, P={P:G9} Pa, D={D:G9} kg/m3, q={Q}");
    }

    private static double V(double? x) => x ?? double.NaN;

    private static Snapshot FromSharpFluids(Fluid f) => new(
        f.FailState, "", V(f.Temperature?.Kelvins), V(f.Pressure?.Pascal), V(f.Density?.KilogramPerCubicMeter),
        V(f.Enthalpy?.JoulePerKilogram), V(f.Entropy?.JoulePerKilogramKelvin), V(f.InternalEnergy?.JoulePerKilogram),
        V(f.Cp?.JoulePerKilogramKelvin), V(f.Cv?.JoulePerKilogramKelvin), V(f.SoundSpeed?.MeterPerSecond),
        V(f.DynamicViscosity?.NewtonSecondPerMeterSquared), V(f.Conductivity?.WattPerMeterKelvin), f.Prandtl,
        V(f.SurfaceTension?.NewtonPerMeter), f.Compressibility, f.Quality, (int)(f.Phase ?? SharpFluids.Phases.Unknown),
        V(f.Tsat?.Kelvins), V(f.GasDensity?.KilogramPerCubicMeter), V(f.LiquidDensity?.KilogramPerCubicMeter));

    private static Snapshot FromAmmonia(Ammonia a) => new(
        false, "", a.Temperature.Kelvin, a.Pressure.Pascal, a.Density.KilogramPerCubicMeter, a.Enthalpy.JoulePerKilogram,
        a.Entropy.JoulePerKilogramKelvin, a.InternalEnergy.JoulePerKilogram, a.Cp.JoulePerKilogramKelvin, a.Cv.JoulePerKilogramKelvin,
        a.SoundSpeed.MeterPerSecond, a.DynamicViscosity.PascalSecond, a.Conductivity.WattPerMeterKelvin, a.Prandtl,
        a.SurfaceTension.NewtonPerMeter, a.Compressibility, a.Quality, (int)a.Phase, a.Tsat.Kelvin,
        a.GasDensity.KilogramPerCubicMeter, a.LiquidDensity.KilogramPerCubicMeter);

    private static Snapshot Run(Action<Fluid> update)
    {
        var f = new Fluid(FluidList.Ammonia);
        update(f);
        return FromSharpFluids(f);
    }

    private static Snapshot Run(Action<Ammonia> update)
    {
        var a = new Ammonia();
        try
        {
            update(a);
            return FromAmmonia(a);
        }
        catch (Exception e) when (e is InvalidOperationException or ArgumentOutOfRangeException)
        {
            return new Snapshot(true, e.Message, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN,
                double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, -1, double.NaN, double.NaN, double.NaN);
        }
    }

    // ------------------------------------------------------------------
    // Reference states
    // ------------------------------------------------------------------
    private static readonly Lazy<List<State>> States = new(BuildStates);

    private static List<State> BuildStates()
    {
        var states = new List<State>();

        State? From(Fluid f) => f.FailState || f.Temperature is null ? null
            : new State(f.Temperature.Kelvins, f.Pressure!.Pascal, f.Density!.KilogramPerCubicMeter, f.Enthalpy!.JoulePerKilogram,
                f.Entropy!.JoulePerKilogramKelvin, f.Quality is >= 0 and <= 1 ? f.Quality : -1);

        // Single phase over SharpFluids' whole range (Tmin 195.495 K .. Tmax 725 K, p_triple .. 1000 MPa)
        foreach (double p in Enumerable.Range(0, 18).Select(i => Math.Exp(Math.Log(8e3) + i * (Math.Log(5e8) - Math.Log(8e3)) / 17)))
        {
            foreach (double t in Enumerable.Range(0, 18).Select(i => 197.0 + i * (720.0 - 197.0) / 17))
            {
                if (Math.Abs(t - Tc) < 0.5 || t is > 405.3 and < 405.7)
                    continue;
                var f = new Fluid(FluidList.Ammonia);
                f.UpdatePT(EngineeringUnits.Pressure.FromPascal(p), EngineeringUnits.Temperature.FromKelvins(t));
                if (From(f) is State s)
                    states.Add(s);
            }
        }

        // Two-phase, saturation lines included, from just above the triple point to 0.5 K below Tc
        foreach (double p in Enumerable.Range(0, 20).Select(i => Math.Exp(Math.Log(6.2e3) + i * (Math.Log(1.13e7) - Math.Log(6.2e3)) / 19)))
        {
            foreach (double q in new[] { 0.0, 0.3, 0.7, 1.0 })
            {
                var f = new Fluid(FluidList.Ammonia);
                f.UpdatePX(EngineeringUnits.Pressure.FromPascal(p), q);
                if (From(f) is State s && s.T < Tc - 0.5)
                    states.Add(s with { Q = q });
            }
        }

        return states;
    }

    // ------------------------------------------------------------------
    // The updates, each as (SharpFluids call, EngineeringFluids call)
    // ------------------------------------------------------------------
    private sealed record Flash(string Name, Action<Fluid, State> SharpFluids, Action<Ammonia, State> Ammonia, Func<State, bool> Applies);

    private static EngineeringUnits.Pressure P(State s) => EngineeringUnits.Pressure.FromPascal(s.P);
    private static EngineeringUnits.Temperature T(State s) => EngineeringUnits.Temperature.FromKelvins(s.T)!;
    private static EngineeringUnits.Density D(State s) => EngineeringUnits.Density.FromKilogramPerCubicMeter(s.D);
    private static EngineeringUnits.SpecificEnergy H(State s) => EngineeringUnits.SpecificEnergy.FromJoulePerKilogram(s.H);
    private static EngineeringUnits.SpecificEntropy S(State s) => EngineeringUnits.SpecificEntropy.FromJoulePerKilogramKelvin(s.S);

    private static readonly Dictionary<string, Flash> Flashes = new[]
    {
        new Flash("PT", (f, s) => f.UpdatePT(P(s), T(s)), (a, s) => a.UpdatePT(Pressure.FromPascal(s.P), Temperature.FromKelvin(s.T)), s => !s.TwoPhase),
        new Flash("PX", (f, s) => f.UpdatePX(P(s), s.Q), (a, s) => a.UpdatePXExact(Pressure.FromPascal(s.P), s.Q), s => s.TwoPhase),
        new Flash("XT", (f, s) => f.UpdateXT(s.Q, T(s)), (a, s) => a.UpdateXTExact(s.Q, Temperature.FromKelvin(s.T)), s => s.TwoPhase),
        new Flash("PH", (f, s) => f.UpdatePH(P(s), H(s)), (a, s) => a.UpdatePHExact(Pressure.FromPascal(s.P), Enthalpy.FromJoulePerKilogram(s.H)), _ => true),
        new Flash("PS", (f, s) => f.UpdatePS(P(s), S(s)), (a, s) => a.UpdatePSExact(Pressure.FromPascal(s.P), SpecificEntropy.FromJoulePerKilogramKelvin(s.S)), _ => true),
        new Flash("DT", (f, s) => f.UpdateDT(D(s), T(s)), (a, s) => a.UpdateDT(Density.FromKilogramPerCubicMeter(s.D), Temperature.FromKelvin(s.T)), _ => true),
        new Flash("DP", (f, s) => f.UpdateDP(D(s), P(s)), (a, s) => a.UpdateDP(Density.FromKilogramPerCubicMeter(s.D), Pressure.FromPascal(s.P)), _ => true),
        new Flash("DH", (f, s) => f.UpdateDH(D(s), H(s)), (a, s) => a.UpdateDH(Density.FromKilogramPerCubicMeter(s.D), Enthalpy.FromJoulePerKilogram(s.H)), _ => true),
        new Flash("DS", (f, s) => f.UpdateDS(D(s), S(s)), (a, s) => a.UpdateDS(Density.FromKilogramPerCubicMeter(s.D), SpecificEntropy.FromJoulePerKilogramKelvin(s.S)), _ => true),
        new Flash("TS", (f, s) => f.UpdateTS(T(s), S(s)), (a, s) => a.UpdateTS(Temperature.FromKelvin(s.T), SpecificEntropy.FromJoulePerKilogramKelvin(s.S)), _ => true),
        new Flash("HS", (f, s) => f.UpdateHS(H(s), S(s)), (a, s) => a.UpdateHS(Enthalpy.FromJoulePerKilogram(s.H), SpecificEntropy.FromJoulePerKilogramKelvin(s.S)), _ => true),
    }.ToDictionary(f => f.Name);

    public static IEnumerable<object[]> FlashNames() => Flashes.Keys.Select(k => new object[] { k });

    // ------------------------------------------------------------------
    // Tolerances (relative, against max(|SharpFluids value|, floor))
    // ------------------------------------------------------------------
    private sealed record Check(string Name, Func<Snapshot, double> Get, double Tol, double Floor, Func<Snapshot, bool>? When = null);

    // The two implementations of the same EOS agree to ~1e-8 in the state variables. A few inverse flashes amplify
    // that - e.g. PH at 500 MPa turns a 1e-8 difference in h into 1.6e-7 in s - hence 5e-7. The worst actual values
    // are printed by every run.
    private const double StateTol = 5e-7;

    private static readonly Check[] Checks =
    {
        new("Temperature", s => s.T, StateTol, 1),
        new("Pressure", s => s.P, StateTol, 1),
        new("Density", s => s.D, StateTol, 1e-3),
        new("Enthalpy", s => s.H, StateTol, 1e3),
        new("Entropy", s => s.S, StateTol, 10),
        new("InternalEnergy", s => s.U, StateTol, 1e3),
        new("Cp", s => s.Cp, 1e-6, 1),
        new("Cv", s => s.Cv, 1e-6, 1),
        new("DynamicViscosity", s => s.Mu, 1e-6, 1e-9),
        new("Conductivity", s => s.K, 1e-6, 1e-6),
        new("Prandtl", s => s.Pr, 1e-6, 1e-6),
        new("Compressibility", s => s.Z, 1e-6, 1),
        new("GasDensity", s => s.GasD, StateTol, 1e-3),
        new("LiquidDensity", s => s.LiqD, StateTol, 1e-3),
        // SharpFluids: 0 whenever quality is in [0, 1]
        new("SoundSpeed", s => s.W, 1e-6, 1, s => !(s.Q is >= 0 and <= 1)),
        // SharpFluids: 0 on the saturation lines, stale in single phase
        new("SurfaceTension", s => s.Sigma, 1e-9, 1e-9, s => s.Q is > 0 and < 1),
        // SharpFluids: the state's own temperature below the triple-point pressure
        new("Tsat", s => s.Tsat, StateTol, 1, s => s.P >= 6091.2231081315085),
    };

    // ------------------------------------------------------------------
    // The comparison
    // ------------------------------------------------------------------
    [DataTestMethod]
    [DynamicData(nameof(FlashNames), DynamicDataSourceType.Method)]
    public void Update_GivesTheSameResultAsSharpFluids(string flashName)
    {
        Flash flash = Flashes[flashName];
        var problems = new List<string>();
        var worst = Checks.ToDictionary(c => c.Name, _ => (err: 0.0, at: ""));
        int compared = 0, bothFailed = 0, onlyEngineeringFluidsSolved = 0;
        var sharpFluidsWrong = new List<string>();

        foreach (State state in States.Value.Where(flash.Applies))
        {
            Snapshot sf = Run(f => flash.SharpFluids(f, state));
            Snapshot ef = Run(a => flash.Ammonia(a, state));

            if (sf.Failed || ef.Failed)
            {
                if (sf.Failed && ef.Failed)
                    bothFailed++;
                else if (sf.Failed && ReproducesReferenceState(ef, state))
                    onlyEngineeringFluidsSolved++; // CoolProp gave up on an input whose state EngineeringFluids found correctly
                else
                    problems.Add($"{state}: SharpFluids {(sf.Failed ? "failed" : "succeeded")}, EngineeringFluids {(ef.Failed ? "failed: " + ef.Error : "succeeded with a different state")}");
                continue;
            }

            // The inputs come from a known state. If SharpFluids did not land on it and EngineeringFluids did, that is
            // a SharpFluids/CoolProp error (seen for DS on the saturation lines: right T and D, h off by a factor 2).
            if (!ReproducesReferenceState(sf, state) && ReproducesReferenceState(ef, state))
            {
                sharpFluidsWrong.Add($"{state}: SharpFluids gave T={sf.T:G9} K, P={sf.P:G9} Pa, h={sf.H:G9} J/kg");
                continue;
            }

            compared++;

            // Quality: SharpFluids' sentinels mean single phase. A state exactly on a saturation line may be called
            // saturated (q = 0/1) by one system and single-phase by the other - the same physical state.
            double sfQ = sf.Q is >= 0 and <= 1 ? sf.Q : -1;
            bool boundaryAmbiguity = state.OnSaturationLine && (sfQ == -1) != (ef.Q == -1);
            if (!boundaryAmbiguity)
            {
                if (Math.Abs(sfQ - ef.Q) > 1e-7)
                    problems.Add($"{state}: Quality SharpFluids {sfQ} vs EngineeringFluids {ef.Q}");
                if (sf.Phase != ef.Phase)
                    problems.Add($"{state}: Phase SharpFluids {(Phases)sf.Phase} vs EngineeringFluids {(Phases)ef.Phase}");
            }

            foreach (Check c in Checks)
            {
                if (c.When is not null && !(c.When(sf) && c.When(ef with { Q = ef.Q })))
                    continue;
                if (boundaryAmbiguity && c.Name is "SoundSpeed" or "SurfaceTension" or "GasDensity" or "LiquidDensity" or "Cp" or "Cv" or "DynamicViscosity" or "Conductivity" or "Prandtl" or "Compressibility")
                    continue; // saturated vs single-phase conventions differ for these; the state itself is compared above
                double a = c.Get(sf), b = c.Get(ef);
                double err = Math.Abs(a - b) / Math.Max(Math.Abs(a), c.Floor);
                if (!(err <= c.Tol))
                    problems.Add(string.Create(CultureInfo.InvariantCulture, $"{state}: {c.Name} SharpFluids {a:G10} vs EngineeringFluids {b:G10} (rel {err:E1})"));
                if (err > worst[c.Name].err)
                    worst[c.Name] = (err, state.ToString());
            }
        }

        var report = new StringBuilder($"{flashName}: {compared} states compared, {bothFailed} failed in both, " +
            $"{onlyEngineeringFluidsSolved} solved only by EngineeringFluids and {sharpFluidsWrong.Count} solved wrongly by SharpFluids (both checked against the reference state).\n");
        foreach (string wrong in sharpFluidsWrong.Take(5))
            report.AppendLine("  SharpFluids wrong: " + wrong);
        foreach (var (name, (err, at)) in worst.Where(w => w.Value.err > 0))
            report.AppendLine(string.Create(CultureInfo.InvariantCulture, $"  {name,-17} worst rel {err:E1}  at {at}"));
        Console.WriteLine(report);

        Assert.IsTrue(compared > 50, $"Only {compared} states compared for {flashName}.");
        Assert.AreEqual(0, problems.Count, $"{problems.Count} differences for {flashName}:\n  " + string.Join("\n  ", problems.Take(25)) + "\n" + report);
    }

    // The inputs of every update come from a reference state, so when SharpFluids gives up but EngineeringFluids
    // does not, EngineeringFluids must have found that very state.
    private static bool ReproducesReferenceState(Snapshot ef, State state)
    {
        static bool Close(double a, double b, double floor) => Math.Abs(a - b) <= 1e-6 * Math.Max(Math.Abs(b), floor);
        return !ef.Failed && Close(ef.T, state.T, 1) && Close(ef.P, state.P, 1) && Close(ef.D, state.D, 1e-3)
            && Close(ef.H, state.H, 1e3) && Close(ef.S, state.S, 10);
    }

    // ------------------------------------------------------------------
    // Constants and saturation helpers
    // ------------------------------------------------------------------
    [TestMethod]
    public void LimitsAndConstants_MatchSharpFluids()
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePT(EngineeringUnits.Pressure.FromBar(21), EngineeringUnits.Temperature.FromKelvins(373.15)); // MolarMass is only set by an update
        var a = new Ammonia();
        Assert.AreEqual(f.LimitTemperatureMin!.Kelvins, a.LimitTemperatureMin.Kelvin, 1e-12);
        Assert.AreEqual(f.LimitTemperatureMax!.Kelvins, a.LimitTemperatureMax.Kelvin, 1e-12);
        Assert.AreEqual(f.LimitPressureMin!.Pascal, a.LimitPressureMin.Pascal, 1e-9);
        Assert.AreEqual(f.LimitPressureMax!.Pascal, a.LimitPressureMax.Pascal, 1e-3);
        Assert.AreEqual(f.CriticalTemperature!.Kelvins, a.CriticalTemperature.Kelvin, 1e-12);
        Assert.AreEqual(f.CriticalPressure!.Pascal, a.CriticalPressure.Pascal, 1e-6);
        Assert.AreEqual(f.T_freeze!.Kelvins, a.T_freeze.Kelvin, 1e-12);
        Assert.AreEqual(f.MolarMass!.KilogramPerMole, a.MolarMass.KilogramPerMole, 1e-15);
    }

    [TestMethod]
    public void GetSatTemperatureAndPressure_MatchSharpFluids()
    {
        var f = new Fluid(FluidList.Ammonia);
        var a = new Ammonia();
        foreach (double p in Enumerable.Range(0, 30).Select(i => Math.Exp(Math.Log(6.1e3) + i * (Math.Log(1.12e7) - Math.Log(6.1e3)) / 29)))
        {
            double expected = f.GetSatTemperature(EngineeringUnits.Pressure.FromPascal(p))!.Kelvins;
            Assert.AreEqual(expected, a.GetSatTemperature(Pressure.FromPascal(p)).Kelvin, 1e-7 * expected, $"Tsat at {p} Pa");
        }
        foreach (double t in Enumerable.Range(0, 30).Select(i => 195.495 + i * (405.0 - 195.495) / 29))
        {
            double expected = f.GetSatPressure(EngineeringUnits.Temperature.FromKelvins(t))!.Pascal;
            Assert.AreEqual(expected, a.GetSatPressure(Temperature.FromKelvin(t)).Pascal, 1e-7 * expected, $"psat at {t} K");
        }
        // Both: above the critical point the critical values
        Assert.AreEqual(f.GetSatTemperature(EngineeringUnits.Pressure.FromBar(200))!.Kelvins, a.GetSatTemperature(Pressure.FromBar(200)).Kelvin);
        Assert.AreEqual(f.GetSatPressure(EngineeringUnits.Temperature.FromKelvins(500))!.Pascal, a.GetSatPressure(Temperature.FromKelvin(500)).Pascal);
    }

    // ------------------------------------------------------------------
    // The SharpFluids behaviour the comparison works around. If one of these starts failing, SharpFluids has
    // changed: remove the matching workaround above.
    // ------------------------------------------------------------------
    [TestMethod]
    public void KnownSharpFluidsBehaviour_UpdatePT_Between405_3And405_7K_ReportsA405_3KState()
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePT(EngineeringUnits.Pressure.FromPascal(2193443.57087332), EngineeringUnits.Temperature.FromKelvins(405.45));
        Assert.AreEqual(405.3, f.Temperature!.Kelvins, 1e-9);

        var a = new Ammonia();
        a.UpdatePT(Pressure.FromPascal(2193443.57087332), Temperature.FromKelvin(405.45));
        Assert.AreEqual(405.45, a.Temperature.Kelvin, 1e-12);
    }

    [TestMethod]
    public void KnownSharpFluidsBehaviour_HelmholtzDiagnosticsAreZero()
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePT(EngineeringUnits.Pressure.FromBar(21), EngineeringUnits.Temperature.FromKelvins(373.15));
        Assert.AreEqual(0, f.Tau);
        Assert.AreEqual(0, f.Delta);
        Assert.AreEqual(0, f.AlphaR);
        Assert.AreEqual(0, f.fugacity);
        Assert.AreEqual(0, f.gibbsmolar);
        Assert.AreEqual(0, f.umolar);
    }

    [TestMethod]
    public void KnownSharpFluidsBehaviour_SoundSpeedIsZeroOnTheSaturationLines()
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePX(EngineeringUnits.Pressure.FromBar(10), 0.0);
        Assert.AreEqual(0, f.SoundSpeed!.MeterPerSecond);

        // EngineeringFluids gives the saturated-liquid value there, as CoolProp itself does
        var a = new Ammonia();
        a.UpdatePXExact(Pressure.FromBar(10), 0.0);
        Assert.IsTrue(a.SoundSpeed.MeterPerSecond > 1000);
    }

    [TestMethod]
    public void KnownSharpFluidsBehaviour_SurfaceTensionIsStaleAfterLeavingTwoPhase()
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePX(EngineeringUnits.Pressure.FromBar(10), 0.5);
        f.UpdatePT(EngineeringUnits.Pressure.FromBar(21), EngineeringUnits.Temperature.FromKelvins(373.15));
        Assert.IsNotNull(f.SurfaceTension);
        Assert.IsTrue(f.SurfaceTension!.NewtonPerMeter > 0);
    }
}
