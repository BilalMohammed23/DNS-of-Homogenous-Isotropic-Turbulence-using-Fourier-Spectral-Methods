# Homogeneous Isotropic Turbulence using Fourier Spectral Methods

This project simulates incompressible homogeneous isotropic turbulence in a three-dimensional $2\pi$-periodic domain using Fourier spectral methods in MATLAB. Homogeneous turbulence assumes spatially uniform statistical properties, while isotropic turbulence assumes direction-independent statistics with equal velocity fluctuation intensities.

The solver is developed for periodic incompressible flow and solves the Fourier-transformed Navier–Stokes equations using spectral differentiation, pressure projection, viscous diffusion, and pseudo-spectral nonlinear term evaluation. De-aliasing is implemented to reduce aliasing errors arising from nonlinear convective products.

The main purpose of this project is to build a MATLAB-based spectral turbulence solver capable of initializing divergence-free velocity fields and analyzing turbulence evolution through quantities such as velocity fields, kinetic energy decay, velocity-derivative skewness, and three-dimensional energy spectra.

## Governing Equations

The solver is based on the incompressible Navier–Stokes equations:

```math
\nabla \cdot \mathbf{u} = 0

```
## Code Structure

```text
main.m
Main driver script for defining the grid, simulation parameters, time loop, diagnostics, and visualization.

uvw_hatk_ic_correct_function.m
Generates divergence-free isotropic turbulence initial conditions in Fourier space.

main_function_code.m
Evaluates the Fourier-space Navier–Stokes right-hand side, including nonlinear convection, pressure gradient, and viscous diffusion.

pressure_poisson.m
Computes the pressure field in Fourier space using the projection/Poisson formulation.

RK4_PsuedoSpectral_Isotropic.m
Advances the velocity Fourier coefficients using the fourth-order Runge–Kutta method.

pre_product_psuedo_spect.m
Applies the 3/2-rule zero-padding before pseudo-spectral nonlinear product evaluation.

deAliased_product_psuedo_spect.m
Computes nonlinear products on the padded grid and truncates them back to the original Fourier grid.

spectrum3D_radial.m
Computes the three-dimensional radial energy spectrum from the Fourier velocity coefficients.

dt.m
Computes the stable time step based on grid spacing, maximum velocity, viscosity, and RK4 stability constant.

```text
