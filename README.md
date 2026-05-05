# Homogenous-Isotropic-Turbulence-using-Fourier-Spectral-Methods
This project simulates incompressible homogeneous isotropic turbulence in a three-dimensional \(2\pi\)-periodic domain using Fourier spectral methods in MATLAB. Homogeneous turbulence assumes spatially uniform statistical properties, while isotropic turbulence assumes direction-independent statistics with equal velocity fluctuation intensities.

The solver is developed for periodic incompressible flow and solves the Fourier-transformed Navier–Stokes equations using spectral differentiation, pressure projection, viscous diffusion, and pseudo-spectral nonlinear term evaluation. De-aliasing is implemented to reduce aliasing errors arising from nonlinear convective products.

The main purpose of this project is to build a MATLAB-based spectral turbulence solver capable of initializing divergence-free velocity fields and analyzing turbulence evolution through quantities such as velocity fields, kinetic energy decay, velocity-derivative skewness, and three-dimensional energy spectra.
