# Homogeneous Isotropic Turbulence using Fourier Spectral Methods

This project simulates incompressible homogeneous isotropic turbulence in a three-dimensional $2\pi$-periodic domain using Fourier spectral methods in MATLAB. Homogeneous turbulence assumes spatially uniform statistical properties, while isotropic turbulence assumes direction-independent statistics with equal velocity fluctuation intensities.

The solver is developed for periodic incompressible flow and solves the Fourier-transformed Navier–Stokes equations using spectral differentiation, pressure projection, viscous diffusion, and pseudo-spectral nonlinear term evaluation. De-aliasing is implemented to reduce aliasing errors arising from nonlinear convective products.

The main purpose of this project is to build a MATLAB-based spectral turbulence solver capable of initializing divergence-free velocity fields and analyzing turbulence evolution through quantities such as velocity fields, kinetic energy decay, velocity-derivative skewness, and three-dimensional energy spectra.

# Numerical Solution of Homogeneous Isotropic Turbulence (DNS — Pseudospectral Method)

---

## Table of Contents

1. [Homogeneous Turbulence](#1-homogeneous-turbulence)
2. [Isotropic Turbulence](#2-isotropic-turbulence)
3. [Reynolds Stress Tensor](#3-reynolds-stress-tensor)
4. [Turbulence Kinetic Energy & Decay](#4-turbulence-kinetic-energy--decay)
5. [Governing Equations — Incompressible Navier–Stokes](#5-governing-equations--incompressible-navierstokes)
6. [Fourier-Space (Spectral) Formulation](#6-fourier-space-spectral-formulation)
7. [Pressure Elimination via Continuity](#7-pressure-elimination-via-continuity)
8. [Pseudospectral Time-Advancement Algorithm](#8-pseudospectral-time-advancement-algorithm)
9. [Initial Conditions — Rogallo Method](#9-initial-conditions--rogallo-method)
10. [Divergence-Free Basis Construction](#10-divergence-free-basis-construction)
11. [Choosing α(k) and β(k) — Isotropy Condition](#11-choosing-αk-and-βk--isotropy-condition)
12. [Energy Spectrum Specification](#12-energy-spectrum-specification)
13. [Taylor Microscale](#13-taylor-microscale)
14. [Simulation Setup — Step-by-Step](#14-simulation-setup--step-by-step)
15. [Key Diagnostics & Results](#15-key-diagnostics--results)

---

## 1. Homogeneous Turbulence

Turbulence is **homogeneous** when its statistics are **independent of spatial position** (translation invariant).

| Configuration | Homogeneous in | Varies in | Numerical treatment |
|---|---|---|---|
| Flat-plate boundary layer | `z` only | `x`, `y` | Finite diff. in `x`, `y` |
| Channel flow | `x`, `z` | `y` (wall-normal) | Periodic in `x`,`z` — finite diff. in `y` |
| Isotropic turbulence | `x`, `y`, `z` | — | Periodic in all directions |

> **Note:** In channel flow, `y` has varying flow properties from wall to center — it is *not* homogeneous in that direction.

---

## 2. Isotropic Turbulence

Isotropic turbulence is a **special case of homogeneous turbulence** with the additional condition that statistics are also **independent of orientation**.

**Implications:**

- Being homogeneous ⟹ Fourier decomposition is valid in all 3 directions
- Statistics do not change in space ⟹ a single **volume-averaged** value characterises the entire domain
- **No mean shear stress** is present (unlike shear-driven flows)
- **Periodic boundary conditions** apply in all three directions `(x, y, z)`
- Velocity values are the same across periodic boundaries

---

## 3. Reynolds Stress Tensor

### Definition

$$R_{ij} = \overline{u_i' u_j'} = R_{ij}(t)$$

where $\overline{(\cdot)}$ denotes a **volume average**, and $u_i'$, $u_j'$ are velocity fluctuations in the $i$ and $j$ directions.

> $R_{ij}$ decays to zero in time (in the DNS sense). In grid-turbulence experiments (e.g., Comte-Bellot & Corrsin, JFM 1971), $R_{ij}$ is a **time average** and decays spatially in $x$.

### Isotropy Conditions

**Diagonal terms** (equal for isotropic turbulence):

$$\overline{u'^2} = \overline{v'^2} = \overline{w'^2}$$

> These are statistically uncorrelated in different directions.

**Off-diagonal terms** (vanish for isotropic turbulence):

$$\overline{u'v'} = \overline{u'w'} = \overline{v'w'} = 0$$

> For shear flows, $\overline{u'v'} \neq 0$.

**Mean velocity** is uniform and commonly set to **zero**.

---

## 4. Turbulence Kinetic Energy & Decay

### Notation

$$q^2 = R_{ii} = \overline{u'^2} + \overline{v'^2} + \overline{w'^2} = 2k$$

where $k$ is the turbulence kinetic energy. Hence $q^2$ is **twice** the TKE.

- $q_0^2$ = initial energy
- $q^2 \to 0$ as $t \to \infty$ (free decay)
- Decay is **slower** at higher Reynolds numbers

### Decay Behaviour

Plotting $q^2 / q_0^2$ vs. time yields a monotonically decreasing curve. The Taylor Reynolds number $R_\lambda$ tracks how quickly this decay proceeds.

---

## 5. Governing Equations — Incompressible Navier–Stokes

### Momentum Equation

$$\frac{\partial u_i}{\partial t} + \frac{\partial (u_i u_j)}{\partial x_j} = -\frac{1}{\rho} \frac{\partial p}{\partial x_i} + \nu \frac{\partial^2 u_i}{\partial x_j \partial x_j}$$

### Continuity (Incompressibility Constraint)

$$\frac{\partial u_i}{\partial x_i} = 0$$

### Fourier Representation (3D Periodic Domain)

Since all three directions are periodic, the velocity field is expanded as a Fourier series:

$$u_i(x, y, z, t) = \sum_{k_x = -N_x/2}^{N_x/2 - 1} \sum_{k_y = -N_y/2}^{N_y/2 - 1} \sum_{k_z = -N_z/2}^{N_z/2 - 1} \hat{u}_i(k_x, k_y, k_z, t) \, e^{i(k_x x + k_y y + k_z z)}$$

where $\vec{k} \cdot \vec{x} = k_i x_i$.

---

## 6. Fourier-Space (Spectral) Formulation

Substituting the Fourier expansion into the N–S equations and invoking **orthogonality** (i.e., Fourier transforming):

$$\frac{\partial \hat{u}_i}{\partial t} + ik_j \widehat{u_i u_j} = -ik_i \hat{P} - \nu(k_1^2 + k_2^2 + k_3^2)\hat{u}_i$$

which simplifies to:

$$\boxed{\frac{\partial \hat{u}_i}{\partial t} = -ik_j \widehat{u_i u_j} - ik_i \hat{P} - \nu k^2 \hat{u}_i} \tag{1}$$

where $k^2 = k_j k_j = |\vec{k}|^2$.

### Continuity in Fourier Space

$$\frac{\partial u_i}{\partial x_i} = 0 \quad \Rightarrow \quad ik_i \hat{u}_i = 0 \quad \Rightarrow \quad k_i \hat{u}_i = 0 \quad \Rightarrow \quad \vec{k} \cdot \hat{\vec{u}} = 0$$

---

## 7. Pressure Elimination via Continuity

Pressure can be **eliminated analytically** by contracting equation (1) with $k_i$:

$$k_i \frac{\partial \hat{u}_i}{\partial t} = -ik_i k_j \widehat{u_i u_j} - ik_i k_j \hat{P} - \nu k^2 \underbrace{k_i \hat{u}_i}_{= 0}$$

Since $\partial(k_i \hat{u}_i)/\partial t = 0$ (continuity), the pressure in Fourier space is:

$$\boxed{\hat{P} = -\frac{k_i k_j}{k^2} \widehat{u_i u_j}}$$

### Substituting $\hat{P}$ back into (1)

$$\frac{\partial \hat{u}_i}{\partial t} = -ik_j \widehat{u_i u_j} + ik_i \frac{k_\ell k_m}{k^2} \hat{u}_\ell \hat{u}_m - \nu k^2 \hat{u}_i \tag{1}$$

where $\ell, m : 1 \to 3$, and the nonlinear term expands as:

$$k_\ell k_m \hat{u}_\ell \hat{u}_m = k_1 k_1 \hat{u}_1 \hat{u}_1 + k_1 k_2 \hat{u}_1 \hat{u}_2 + k_1 k_3 \hat{u}_1 \hat{u}_3 + k_2 k_1 \hat{u}_2 \hat{u}_1 + \cdots$$

> **Practical note:** It is cleaner to carry $\hat{P}$ as a **separate variable** and solve for it independently — this is also required if pressure contours need to be plotted.

---

## 8. Pseudospectral Time-Advancement Algorithm

Given initial conditions, each time step proceeds as:

```
1. Fourier transform  uᵢ → ûᵢ          (u, v, w)
2. Compute uᵢuⱼ in physical space, then Fourier transform with de-aliasing   [Pseudo-spectral]
3. Integrate the ODE for ûᵢ forward in time
4. Inverse Fourier transform  ûᵢ → uᵢ
```

### Divergence-Free Check

At **every time step**, verify:

$$k_i \hat{u}_i = 0$$

This confirms the solution remains divergence-free throughout the simulation.

---

## 9. Initial Conditions — Rogallo Method

We need to generate $\hat{u}_i(\vec{k})$ at $t = 0$. Three conditions must be satisfied:

| # | Requirement | Description |
|---|---|---|
| i | **Divergence-free** | $\vec{k} \cdot \hat{\vec{u}} = 0$ (continuity) |
| ii | **Real velocity field** | $\hat{u}_i(\vec{k}) = \hat{u}_i^*(-\vec{k})$ (Hermitian symmetry) |
| iii | **Isotropic** (statistically) | Statistics independent of orientation |

> This approach is due to **Rogallo** — the classical method for generating turbulence initial conditions. At the end, we need the Fourier coefficients of the velocity field.

**Satisfying condition ii (Reality):**
Generate $\hat{u}_i$ in a **half-plane** and ensure Hermitian symmetry. $\hat{\vec{u}}$ can be represented as components of $e_1'$ and $e_2'$:

$$\hat{\vec{u}} = \alpha(k) e_1' + \beta(k) e_2' \quad \text{where } \vec{k} = k_1 e_1 + k_2 e_2 + k_3 e_3$$

Here $\alpha(k)$, $\beta(k)$ are the magnitudes to be determined from the energy spectrum (depend on spectrum).

---

## 10. Divergence-Free Basis Construction

Since $\vec{k} \cdot \hat{\vec{u}} = 0$, we require $\vec{k} \perp \hat{\vec{u}}$.

**Strategy:** Generate Fourier coefficients in the plane $\perp \vec{k}$ — this **guarantees** divergence-freedom.

### New Orthonormal Basis $\{e_1', e_2', e_3'\}$

Define a rotated basis where $e_3' \parallel \vec{k}$:

$$k_1 e_1 + k_2 e_2 + k_3 e_3 = k \cdot e_3'$$

Assume $e_1' \perp e_3$ (i.e., $e_1' \cdot e_3 = 0$, meaning $e_1' \perp e_3'$). Using orthogonality:

$$e_1' \cdot e_3' = 0 \quad \Rightarrow \quad Ak_1 + Bk_2 = 0 \quad \Rightarrow \quad B = -A\frac{k_1}{k_2}$$

Solving and normalising (Pythagorean theorem):

$$\boxed{e_1' = \frac{k_2 e_1 - k_1 e_2}{\sqrt{k_1^2 + k_2^2}}}$$

The third basis vector is obtained by cross product:

$$e_2' = e_3' \times e_1'$$

$$k\sqrt{k_1^2 + k_2^2} \; e_2' = k_1 k_3 e_1 + k_2 k_3 e_2 - (k_1^2 + k_2^2) e_3$$

> $e_1', e_2', e_3'$ are a new basis where $e_3' \parallel \vec{k}$ and $k = |\vec{k}| = \sqrt{k_i k_i}$.

### Rogallo Initial Condition (Final Form)

Substituting $e_1'$, $e_2'$ and knowing $\alpha$, $\beta$:

$$\hat{\vec{u}} = \underbrace{\left[\frac{\alpha k k_2 + \beta k_1 k_3}{k\sqrt{k_1^2 + k_2^2}}\right]}_{\hat{u} \; @\; t=0} e_1 \; + \; \underbrace{\left[\frac{\beta k_2 k_3 - \alpha k k_1}{k\sqrt{k_1^2 + k_2^2}}\right]}_{\hat{v} \; @\; t=0} e_2 \; - \; \underbrace{\left[\frac{\beta\sqrt{k_1^2 + k_2^2}}{k}\right]}_{\hat{w} \; @\; t=0} e_3$$

> **Exercise:** Re-derive this expression independently to verify.

---

## 11. Choosing α(k) and β(k) — Isotropy Condition

$\alpha$ and $\beta$ depend only on the **magnitude** $k = |\vec{k}|$. Set:

$$\alpha(k) = \sqrt{\frac{E(k)}{4\pi k^2}} \; e^{i\theta_1} \cos\phi$$

$$\beta(k) = \sqrt{\frac{E(k)}{4\pi k^2}} \; e^{i\theta_2} \sin\phi$$

where $\theta_1$, $\theta_2$, $\phi$ are **random numbers drawn from $[0, 2\pi]$**, and $E(k)$ is the prescribed energy spectrum.

### Verification

This choice ensures:

$$(\alpha\alpha^* + \beta\beta^*) \cdot 4\pi k^2 = E(k)$$

**Energy check** (integrating over all wavenumbers in spherical shells):

$$\int \hat{\vec{u}} \cdot \hat{\vec{u}}^* \, d^3\vec{k} = \int_0^\infty (\alpha\alpha^* + \beta\beta^*) 4\pi k^2 \, dk = \int_0^\infty E(k) \, dk = \text{net kinetic energy}$$

---

## 12. Energy Spectrum Specification

### Special Cases

| Condition | Consequence |
|---|---|
| $k = 0$ ($k_1 = k_2 = k_3 = 0$) | Mean velocity set to **zero** |
| $k_1^2 + k_2^2 = 0$ ($k_1 = k_2 = 0$) | Continuity ⟹ $\hat{u}_3 = 0$ for all $k_3$; set $\hat{u}_1 = \alpha(k)$, $\hat{u}_2 = \beta(k)$ |

### Prescribed Form

One commonly used form:

$$\boxed{E(k) = 16\sqrt{\frac{2}{\pi}} \; \frac{u_0^2}{k_0} \left(\frac{k}{k_0}\right)^4 e^{-2(k/k_0)^2}}$$

| Parameter | Description |
|---|---|
| $k_0$ | Wavenumber at which $E(k)$ is maximum (large scale of the flow) |
| $u_0$ | RMS velocity |

**Normalisation:**

$$\int_0^\infty E(k) \, dk = \frac{3}{2} u_0^2$$

> **Caution:** Higher wavenumbers generally have more numerical error (smaller wavelength → resolution may be insufficient; also prone to aliasing errors).

---

## 13. Taylor Microscale

### General Definition

$$\lambda_\alpha = \sqrt{\frac{\overline{u_\alpha'^2}}{\left(\dfrac{\partial u_\alpha'}{\partial x_\alpha}\right)^2}}$$

This has units of length.

### For the Prescribed Spectrum

$$\lambda = \frac{2}{k_0}$$

For isotropic turbulence: $\lambda_1 = \lambda_2 = \lambda_3$ (statistically). Add all three and divide by 3 to get $\lambda$.

### Taylor Reynolds Number

$$Re_\lambda = \frac{u_0 \lambda}{\nu} = \frac{u_0}{\nu} \cdot \frac{2}{k_0} \quad \Rightarrow \quad \boxed{\nu = \frac{2u_0}{k_0 \cdot Re_\lambda}}$$

### Verification Check

Once the velocity field is generated, compute $\lambda$ using the general definition and compare with $\lambda = 2/k_0$. Agreement confirms the initial condition is correctly set up.

---

## 14. Simulation Setup — Step-by-Step

```
Step 1 — Choose k₀    e.g., k₀ = 5
          (Range of k is [0, N/2] — Nyquist limit)
          If k₀ too small  → statistics become a problem
          If k₀ too large  → resolution becomes a problem

Step 2 — Choose u₀    e.g., u₀ = 1  (non-dimensional)

Step 3 — Choose ν to get the desired Reynolds number
          Re_λ = u₀λ/ν   →   ν = 2u₀ / (k₀ · Re_λ)
          Prescribe Re_λ, then choose u₀, k₀ as 1, 5.

Step 4 — For each (k₁, k₂, k₃), compute k = √(kᵢkᵢ)
          → E(k) determined
          → α(k), β(k) determined
```

---

## 15. Key Diagnostics & Results

### Quantities to Monitor

| Diagnostic | Expression | Notes |
|---|---|---|
| Normalised TKE | $q^2 / q_0^2$ vs. time | Decays monotonically |
| Taylor Reynolds number | $R_\lambda$ vs. time | Tracks decay rate |
| Velocity derivative skewness | $S_u = \overline{u'^3_{,x}} \, / \, [\overline{u'^2_{,x}}]^{3/2}$ | Should settle to $-0.4$ to $-0.6$ for isotropic turbulence |
| Radial energy spectrum | $E(k)$ projected onto spherical shells in $k$-space | Plot at intervals of $\tau_b = 1/(u_0 k_0 / 2)$ |

### Expected Results

- $q^2/q_0^2$ decreases monotonically with time
- Decay is **slower** at higher $Re_\lambda$
- Skewness $S = (S_1 + S_2 + S_3)/3 \approx -0.4$ to $-0.6$

### Reference

> Comte-Bellot, G. & Corrsin, S. (1971). *Simple Eulerian time correlation of full- and narrow-band velocity signals in grid-generated isotropic turbulence.* Journal of Fluid Mechanics.

---

*Notes based on lecture material — DNS of Homogeneous Isotropic Turbulence using the Pseudospectral method (Rogallo, 1981).*


To run the code, Open MATLAB,Keep all files in the same folder, Run main.m

Input parameters are made simple just N and Reynold's number. If you wanna increase the Reynolds ' number, then it's recommended to increase the grid size as well.
