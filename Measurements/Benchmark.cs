using BenchmarkDotNet.Attributes;
using SharpFluids;
using EngineeringFluids.Fluids;
using EngineeringUnits;

namespace Benchmark;

[ShortRunJob]
[MemoryDiagnoser]
public class Benchy
{
    private static readonly Pressure pressure = Pressure.FromBar(21);
    private static readonly Temperature temperature = Temperature.FromDegreeCelsius(100);

    private static readonly double pressured = 2100000; // Pressure.FromBar(21).SI;
    private static readonly double temperatured = 373.15; // Temperature.FromDegreeCelsius(100).SI;

    private static readonly Fluid SharpFluid = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluid = new AmmoniaDouble();
    private static readonly AbstractState CoolPropState = AbstractState.factory("HEOS", "Ammonia");

    [Benchmark]
    public void SharpFluidupdatePT()
    {
        SharpFluid.UpdatePT(pressure, temperature);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePT()
    {
        CoolPropState.update(input_pairs.PT_INPUTS, pressured, temperatured);
    }

    [Benchmark]
    public void EngineeringUpdatePT()
    {
        EngineeringFluid.UpdatePT(pressured, temperatured);
    }

    // ---------------------------------------------------------------------
    // UpdatePX (pressure + vapor-quality flash): same three-way comparison as UpdatePT
    // above, at one representative two-phase point. Unlike UpdatePT, the Fast UpdatePX
    // path has no Newton iteration (it reads Tsat/rhoL/rhoV straight from fixed-cost
    // ancillary polynomials - see UpdatePXFast.cs), so its cost doesn't vary by region
    // the way UpdatePT's does; one point is enough to characterize it.
    // ---------------------------------------------------------------------

    private static readonly Fluid SharpFluidPX = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidPX = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStatePX = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Pressure pxPressure = Pressure.FromPascal(1_000_000);
    private static readonly double pxPressured = 1_000_000;
    private static readonly double pxQuality = 0.5;

    [Benchmark]
    public void SharpFluidUpdatePX()
    {
        SharpFluidPX.UpdatePX(pxPressure, pxQuality);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePX()
    {
        CoolPropStatePX.update(input_pairs.PQ_INPUTS, pxPressured, pxQuality);
    }

    [Benchmark]
    public void EngineeringUpdatePX()
    {
        EngineeringFluidPX.UpdatePX(pxPressured, pxQuality);
    }

    // ---------------------------------------------------------------------
    // UpdatePH (pressure + enthalpy flash): same three-way comparison, at both a single-phase
    // point and a two-phase point. Unlike UpdatePX, the Fast UpdatePH path's cost genuinely
    // varies by case: the two-phase branch is an O(1) ancillary lookup (same class as UpdatePX),
    // but the single-phase branch runs an outer temperature root-find with UpdatePT as the
    // inner evaluator (see UpdatePHFast.cs), so it costs several UpdatePT-and-property-read
    // iterations rather than one lookup.
    // ---------------------------------------------------------------------

    private static readonly Fluid SharpFluidPH = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidPH = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStatePH = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Pressure phPressure = Pressure.FromPascal(2_100_000);
    private static readonly SpecificEnergy phEnthalpy = SpecificEnergy.FromJoulePerKilogram(1_792_554.16228943);
    private static readonly double phPressured = 2_100_000;
    private static readonly double phEnthalpyd = 1_792_554.16228943;

    private static readonly Fluid SharpFluidPH2P = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidPH2P = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStatePH2P = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Pressure ph2pPressure = Pressure.FromPascal(1_000_000);
    private static readonly SpecificEnergy ph2pEnthalpy = SpecificEnergy.FromJoulePerKilogram(1_045_846.5);
    private static readonly double ph2pPressured = 1_000_000;
    private static readonly double ph2pEnthalpyd = 1_045_846.5;

    [Benchmark]
    public void SharpFluidUpdatePH()
    {
        SharpFluidPH.UpdatePH(phPressure, phEnthalpy);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePH()
    {
        CoolPropStatePH.update(input_pairs.HmassP_INPUTS, phEnthalpyd, phPressured);
    }

    [Benchmark]
    public void EngineeringUpdatePH()
    {
        EngineeringFluidPH.UpdatePH(phPressured, phEnthalpyd);
    }

    [Benchmark]
    public void SharpFluidUpdatePH_TwoPhase()
    {
        SharpFluidPH2P.UpdatePH(ph2pPressure, ph2pEnthalpy);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePH_TwoPhase()
    {
        CoolPropStatePH2P.update(input_pairs.HmassP_INPUTS, ph2pEnthalpyd, ph2pPressured);
    }

    [Benchmark]
    public void EngineeringUpdatePH_TwoPhase()
    {
        EngineeringFluidPH2P.UpdatePH(ph2pPressured, ph2pEnthalpyd);
    }

    // ---------------------------------------------------------------------
    // UpdatePS (pressure + entropy flash): same three-way comparison and cost shape as
    // UpdatePH above (two-phase is an O(1) ancillary lookup, single-phase is a warm-started
    // outer temperature root-find - see UpdatePSFast.cs). Uses the SAME physical states as
    // the UpdatePH benchmarks (P=2.1MPa gas / P=1MPa q=0.5 mixture), just entropy instead of
    // enthalpy as the second input, so the two sets of numbers are directly comparable.
    // ---------------------------------------------------------------------

    private static readonly Fluid SharpFluidPS = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidPS = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStatePS = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Pressure psPressure = Pressure.FromPascal(2_100_000);
    private static readonly double psPressured = 2_100_000;
    private static readonly SpecificEntropy psEntropy = EntropyAtPT(2_100_000, 373.15);
    private static readonly double psEntropyd = psEntropy.JoulePerKilogramKelvin;

    private static readonly Fluid SharpFluidPS2P = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidPS2P = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStatePS2P = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Pressure ps2pPressure = Pressure.FromPascal(1_000_000);
    private static readonly double ps2pPressured = 1_000_000;
    private static readonly SpecificEntropy ps2pEntropy = EntropyAtPX(1_000_000, 0.5);
    private static readonly double ps2pEntropyd = ps2pEntropy.JoulePerKilogramKelvin;

    private static SpecificEntropy EntropyAtPT(double pPa, double tK)
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePT(Pressure.FromPascal(pPa), Temperature.FromKelvin(tK));
        return f.Entropy!;
    }

    private static SpecificEntropy EntropyAtPX(double pPa, double quality)
    {
        var f = new Fluid(FluidList.Ammonia);
        f.UpdatePX(Pressure.FromPascal(pPa), quality);
        return f.Entropy!;
    }

    [Benchmark]
    public void SharpFluidUpdatePS()
    {
        SharpFluidPS.UpdatePS(psPressure, psEntropy);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePS()
    {
        CoolPropStatePS.update(input_pairs.PSmass_INPUTS, psPressured, psEntropyd);
    }

    [Benchmark]
    public void EngineeringUpdatePS()
    {
        EngineeringFluidPS.UpdatePS(psPressured, psEntropyd);
    }

    [Benchmark]
    public void SharpFluidUpdatePS_TwoPhase()
    {
        SharpFluidPS2P.UpdatePS(ps2pPressure, ps2pEntropy);
    }

    [Benchmark]
    public void CoolPropDirectUpdatePS_TwoPhase()
    {
        CoolPropStatePS2P.update(input_pairs.PSmass_INPUTS, ps2pPressured, ps2pEntropyd);
    }

    [Benchmark]
    public void EngineeringUpdatePS_TwoPhase()
    {
        EngineeringFluidPS2P.UpdatePS(ps2pPressured, ps2pEntropyd);
    }

    // ---------------------------------------------------------------------
    // UpdateTX (temperature + vapor-quality flash): even cheaper than UpdatePX - T is already
    // the ancillaries' independent variable, so there's no P->T inversion at all (see
    // UpdateTXFast.cs). Uses SharpFluids' own UpdateXT naming (quality, T) - that's the only
    // signature the package exposes.
    // ---------------------------------------------------------------------

    private static readonly Fluid SharpFluidTX = new Fluid(FluidList.Ammonia);
    private static readonly AmmoniaDouble EngineeringFluidTX = new AmmoniaDouble();
    private static readonly AbstractState CoolPropStateTX = AbstractState.factory("HEOS", "Ammonia");

    private static readonly Temperature txTemperature = Temperature.FromKelvin(280.0);
    private static readonly double txTemperatured = 280.0;
    private static readonly double txQuality = 0.5;

    [Benchmark]
    public void SharpFluidUpdateTX()
    {
        SharpFluidTX.UpdateXT(txQuality, txTemperature);
    }

    [Benchmark]
    public void CoolPropDirectUpdateTX()
    {
        CoolPropStateTX.update(input_pairs.QT_INPUTS, txQuality, txTemperatured);
    }

    [Benchmark]
    public void EngineeringUpdateTX()
    {
        EngineeringFluidTX.UpdateTX(txTemperatured, txQuality);
    }

    // ---------------------------------------------------------------------
    // UpdateTXExact: same operation as UpdateTX above, but solving the real phase-equilibrium
    // conditions against the EOS (SolveAtTFast) instead of reading the pre-fitted ancillary
    // curves - see UpdateTXExact.cs. Same physical state as the UpdateTX benchmark, so the two
    // are directly comparable: this is the accuracy/speed tradeoff for the exact version.
    // ---------------------------------------------------------------------

    private static readonly AmmoniaDouble EngineeringFluidTXExact = new AmmoniaDouble();

    [Benchmark]
    public void EngineeringUpdateTXExact()
    {
        EngineeringFluidTXExact.UpdateTXExact(txTemperatured, txQuality);
    }

    // ---------------------------------------------------------------------
    // Region-specific cases for the Fast (AmmoniaDouble) path only.
    //
    // UpdatePT's Newton solver needs a different number of iterations depending
    // on where (T,P) falls relative to the saturation dome and the critical
    // point (measured via a full-domain sweep while tuning the initial-guess
    // logic): ~2-3 iterations for superheated gas far from the dome, ~4-5 near
    // it or in the liquid region, and up to an order of magnitude more right
    // near the critical point. The single point above (superheated gas, far
    // from the dome) only exercises the cheapest case, so a regression that
    // only shows up in liquid/supercritical/near-critical states would be
    // invisible there. These cases give each region its own tracked number.
    // ---------------------------------------------------------------------

    private static readonly AmmoniaDouble GasNearDomeFluid = new AmmoniaDouble();
    private static readonly double gasNearDomePressure = 5_630_000; // ~0.9x Psat(373.15K): gas, close to the dew line
    private static readonly double gasNearDomeTemperature = 373.15;

    private static readonly AmmoniaDouble LiquidFluid = new AmmoniaDouble();
    private static readonly double liquidPressure = 2_100_000;
    private static readonly double liquidTemperature = 280.0; // subcooled liquid

    private static readonly AmmoniaDouble SupercriticalTypicalFluid = new AmmoniaDouble();
    private static readonly double supercriticalTypicalPressure = 15_000_000; // ~1.3x Pc, comfortably away from Tc
    private static readonly double supercriticalTypicalTemperature = 500.0;

    private static readonly AmmoniaDouble SupercriticalNearCriticalFluid = new AmmoniaDouble();
    private static readonly double supercriticalNearCriticalPressure = 14_158_915.68768298; // worst case found while sweeping the supercritical region
    private static readonly double supercriticalNearCriticalTemperature = 406.56; // 1 K above Tc=405.56K

    [Benchmark]
    public void Gas_NearDome_Fast()
    {
        GasNearDomeFluid.UpdatePT(gasNearDomePressure, gasNearDomeTemperature);
    }

    [Benchmark]
    public void Liquid_Fast()
    {
        LiquidFluid.UpdatePT(liquidPressure, liquidTemperature);
    }

    [Benchmark]
    public void Supercritical_Typical_Fast()
    {
        SupercriticalTypicalFluid.UpdatePT(supercriticalTypicalPressure, supercriticalTypicalTemperature);
    }

    [Benchmark]
    public void Supercritical_NearCritical_Fast()
    {
        SupercriticalNearCriticalFluid.UpdatePT(supercriticalNearCriticalPressure, supercriticalNearCriticalTemperature);
    }
}