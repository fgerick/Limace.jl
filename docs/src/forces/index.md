# Forces

Different forces included in `Limace.jl`. 

So far, there is no heat equation in `Limace.jl`, so the forces include only the mass-term (or `Limace.inertial`), the Coriolis force (`Limace.coriolis`), the Lorentz force (`Limace.lorentz`) and diffusion (`Limace.diffusion`). 
In the induction equation we have the magnetic induction term $\nabla\times\mathbf{u}\times\mathbf{B}$ (`Limace.induction`). Some details for each of these forcings is given on the respective pages.