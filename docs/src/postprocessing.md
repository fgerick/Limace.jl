
# Postprocessing

A non-exhaustive collection of spectral to spatial discretization and postprocessing routines.

## Spatial discretization

In the `Limace.Discretization` submodule, some methods to discretize `BasisElement` and eigenvectors are provided.

```@autodocs
Modules = [Limace.Discretization]
```

## Filtering of spectrum

```@docs
Limace.Processing.eigenvalue_filter
Limace.Processing.numerical_filter
Limace.Processing.observability_filter
```

## Energy and spectra
```@docs
Limace.Processing.energies
Limace.Processing.epeak_etrunc_cartesian
Limace.Processing.lmn_n
Limace.Processing.spectrum
Limace.Processing.spectrum_cartesian
```
