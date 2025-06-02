
# Postprocessing

A non-exhaustive collection of spectral to spatial discretization and postprocessing routines.

## Spatial discretization

In the `Limace.Discretization` submodule, some methods to discretize [BasisElement](@ref) and eigenvectors are provided.

A fully manual discretization routine is given by
```@docs
Limace.Discretization.discretize
```
To evaluate one [BasisElement](@ref) at a given point $(r,\theta,\phi)$, we can do
```@example processing
using Limace #hide
using Limace.Discretization: discretize
u = Inviscid(10)
b = BasisElement(u, Poloidal, (2,0,1))
r,θ,ϕ = 0.9, π/4, 3π/2
discretize(b, r,θ,ϕ)
```

The same also works for a vector of [BasisElement](@ref)s.
```@example processing
bs = [BasisElement(u, Poloidal, (2,0,1)), BasisElement(u, Poloidal, (2,0,2))]
discretize(bs, r,θ,ϕ)
```

To discretize an eigenvector using [Limace.Discretization.discretize](@ref), consider the following example of the Malkus modes
```@example processing
N = 6
Le = 1e-1

u = Inviscid(N)
b = PerfectlyConducting(N)

bases = [u,b]

B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15))
forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u, b, B0), Limace.InductionB0(b,u,B0)]

problem = LimaceProblem(bases, forcings)
Limace.assemble!(problem)
Limace.solve!(problem)
nothing #hide
```

We can chose a single eigenvector `x` and define some arbitrary spatial grid for `r`, `θ` and `ϕ`:
```@example processing
x = problem.sol.vectors[:,1]
nr,nθ,nϕ = 20,30,40
r = range(0.01,0.99, length=nr)
θ = range(0.01,π-0.01,length=nθ)
ϕ = range(0,2π, length=nϕ)
ur,uθ,uϕ, br, bθ, bϕ = discretize(x, u, b, r, π/2 .- θ, ϕ)
size(ur)
```

For more performant transforms to spatial space, it is recommended to use [Limace.Discretization.spectospat](@ref), 
which relies on fast vector spherical harmonic transforms provided through [SHTns.jl](https://github.com/fgerick/SHTns.jl).
```@docs
Limace.Discretization.spectospat
```

To discretize the eigenvector `x` on a Gauss-Legendre grid in `r` and `θ`, as well as equidistant grid in `ϕ`, we can call `spectospat`:
```@example processing
using Limace.Discretization: spectospat
ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(x, u, b, nr, nθ, nϕ)
nothing #hide
```

## Filtering of spectrum

When the background state couples different spectral degrees (may be in radius, spherical harmonic degree or order), the solutions are only truncated approximations of the true solution. 
This is not the case, for example, for the [Inviscid inertial modes](@ref "Inviscid inertial modes in the sphere") or the [Malkus modes](@ref "Malkus modes in the sphere"), which are fully resolved at one truncation degree and are therefore perfectly converged.
In all other cases, the computed eigen spectrum will contain some unconverged solutions that need to be filtered out before further analysis.

To check convergence of the solutions, `Limace.jl` provides the following filters 
```@docs
Limace.Processing.eigenvalue_filter
Limace.Processing.eigenvector_filter
Limace.Processing.numerical_filter
```

After numerical convergence is validated, one can further filter the spectrum using for example a filter that determines the observability of a mode. 
In the geomagnetic context, observability of a mode may be determined by the peak spherical harmonic degree of the poloidal magnetic field of a mode, 
which should lie below the maximum resolvable degree of the considered geomagnetic model (e.g. $l\leq 17$ for the secular variation in the CHAOS-7 model).

```@docs
Limace.Processing.observability_filter
```

## Energy and spectra

To compute the (kinetic and magnetic) energies of the modes, one can use
```@docs
Limace.Processing.energies
```

See [the example usage of `energies` in the Malkus mode example](@ref "Spectrum of modes"), which shows how to plot a frequency vs. energy-ratio spectrum.

`Limace.jl` also provides routines to compute the spectrum of an eigenvector as a function of spherical harmonic degree, spherical harmonic order, radial degree, or Cartesian degree.
```@docs
Limace.Processing.spectrum
Limace.Processing.spectrum_cartesian
```

For example, we can compute the spectra of the modes and plot them for one solution like this:
```@example processing
ls, p_up, p_ut, p_bp, p_bt = Limace.Processing.spectrum(problem, lmn=1)

using CairoMakie

let
	f = Figure()
	ax = Axis(f[1,1], xlabel=L"l", ylabel=L"p(l)", yscale=log10)
	for (spec,label) in zip((p_up, p_ut, p_bp, p_bt),(L"u_p", L"u_t", L"b_p", L"b_t"))
		scatterlines!(ax, ls, spec[:,1].+eps(); label)
	end
	axislegend(ax)
	f
end
```

and to do the same as a function of Cartesian polynomial degree

```@example processing

p_up, p_ut, p_bp, p_bt = Limace.Processing.spectrum_cartesian(problem.sol.vectors, problem.bases...)
ns = 1:problem.bases[1].N

let
	f = Figure()
	ax = Axis(f[1,1], xlabel=L"l", ylabel=L"p(l)", yscale=log10)
	for (spec,label) in zip((p_up, p_ut, p_bp, p_bt),(L"u_p", L"u_t", L"b_p", L"b_t"))
		scatterlines!(ax, ns, spec[:,1].+eps(); label)
	end
	axislegend(ax)
	f
end
```

## Misc functions

Some convenient functions that are used for other postprocessing routines.

```@docs
Limace.Processing.epeak_etrunc_cartesian
Limace.Processing.lmn_n
```