# EngineeringFluids

Fast, unit-typed thermodynamic and transport properties for engineering fluids in .NET.

Built on [EngineeringUnits.Fast](https://github.com/MadsKirkFoged/EngineeringUnits): every input and output is a
typed quantity (`Pressure`, `Temperature`, `Enthalpy`, ...), so there are no unit mix-ups and no conversions by hand.

Currently supports **Ammonia (R717)** using the Gao et al. (2020) Helmholtz equation of state, with results that
match [CoolProp](http://www.coolprop.org/) / [SharpFluids](https://github.com/MadsKirkFoged/SharpFluids).

## Install

```
dotnet add package EngineeringFluids
```

## Quick start

```csharp
using EngineeringFluids.Fluids;
using EngineeringUnits.Fast;

var nh3 = new Ammonia();

// Single-phase state from pressure and temperature
nh3.UpdatePT(Pressure.FromBar(10), Temperature.FromDegreeCelsius(50));

Console.WriteLine(nh3.Density.KilogramPerCubicMeter);
Console.WriteLine(nh3.Enthalpy.JoulePerKilogram);
Console.WriteLine(nh3.Cp.JoulePerKilogramKelvin);
Console.WriteLine(nh3.DynamicViscosity.PascalSecond);
Console.WriteLine(nh3.Phase);

// Two-phase state from pressure and vapor quality (0 = saturated liquid, 1 = saturated vapor)
nh3.UpdatePX(Pressure.FromBar(5), 0.3);

Console.WriteLine(nh3.Temperature.DegreeCelsius);   // saturation temperature
Console.WriteLine(nh3.Quality);
```

## State updates

| Method | Inputs |
|---|---|
| `UpdatePT` | pressure, temperature (optionally a phase hint) |
| `UpdatePH` / `UpdatePHExact` | pressure, enthalpy |
| `UpdatePS` / `UpdatePSExact` | pressure, entropy |
| `UpdatePX` / `UpdatePXExact` | pressure, vapor quality |
| `UpdateTX` / `UpdateTXExact` | temperature, vapor quality |
| `UpdateTS` | temperature, entropy |
| `UpdateHS` | enthalpy, entropy |
| `UpdateDT` | density, temperature |
| `UpdateDP` | density, pressure |
| `UpdateDH` / `UpdateDS` | density, enthalpy / entropy |

The plain methods read the saturation dome from fast ancillary fits. The `Exact` variants solve the phase
equilibrium against the equation of state instead: a bit slower, but they agree with CoolProp to about 1e-10 and
work all the way up to the critical point. Use `Exact` when you need CoolProp-identical results.

## Properties

Pressure, temperature, density, enthalpy, entropy, internal energy, Cp, Cv, speed of sound, viscosity, thermal
conductivity, Prandtl number, surface tension, fugacity, compressibility, phase, quality, saturation
temperature/pressure, critical and triple-point data, and validity limits.

All properties of a state are computed from a single cached equation-of-state evaluation, so reading many
properties costs about the same as reading one. Once a state is set, it can be read from multiple threads at the same time.

## License

MIT
