# High-Precision Numerical Orbit Integration &amp; Non-Gravitational Perturbations for Interstellar Comet 3I/ATLAS

Welcome to the official GitHub repository for the high-precision numerical orbit integration and non-gravitational perturbation model of the interstellar hyperbolic comet **3I/ATLAS**.

This repository contains the Python simulation source code (`3iatlas_simulation_code-v2.py`) validating the **3iatlas** platform and accompanying the peer-reviewed research manuscript submitted to the **Icarus (Elsevier)** journal.

---

## 🌌 Overview &amp; Physical Model

Interstellar objects traversing the inner Solar System at extreme velocities exhibit minute non-gravitational acceleration due to solar radiation pressure and outgassing thrust.

Based on observational parameters from NASA/JPL, the heliocentric trajectory of **3I/ATLAS** is integrated using an adaptive 4th/5th-order Runge-Kutta integrator (`RK45` via `scipy.integrate.solve_ivp`) with strict error tolerances (`rtol=1e-9`, `atol=1e-9`).

The physical model includes:

1. **Newtonian Solar Gravitation**: Monolithic solar attractor ($\\mu\_\\odot = 1.327124 \\times 10^{11} \\text{ km}^3/\\text{s}^2$).
2. **Solar Radiation Pressure**: Modeled for a porous, low-density nucleus with:  
  * Nucleus Radius ($R$): $1.3\\text{ km}$ ($1300\\text{ m}$)
  * Bulk Density ($\\rho$): $550\\text{ kg/m}^3$
  * Nucleus Mass ($M\_0$): $\\approx 5.06 \\times 10^{12}\\text{ kg}$
  * Cross-sectional Area ($A$): $\\approx 5.31 \\times 10^6\\text{ m}^2$
3. **Asymmetric Outgassing Thrust**: Active when heliocentric distance $r &lt; 2.5\\text{ AU}$, assuming an $8^\\circ$ jet misalignment relative to the Sun-comet radial vector, producing non-conservative orbital energy shifts and a cumulative outgassing drift of **$236.36\\text{ meters}$**.

---

## 📊 Summary of Physical and Simulation Results

| Parameter                         | Value                                                                          | Description                              |
| --------------------------------- | ------------------------------------------------------------------------------ | ---------------------------------------- |
| **Hyperbolic Eccentricity ($e$)** | $6.1414$                                                                       | Orbital eccentricity (NASA JPL)          |
| **Excess Speed ($v\_\\infty$)**   | $64.04\\text{ km/s}$                                                           | Asymptotic inbound speed relative to Sun |
| **Simulated Perihelion ($q$)**    | $0.9388\\text{ AU}$ ($1.404 \\times 10^8\\text{ km}$)                          | Minimum heliocentric approach distance   |
| **Max Non-Gravitational Accel.**  | $6.48 \\times 10^{-15}\\text{ km/s}^2$ ($6.48 \\times 10^{-12}\\text{ m/s}^2$) | Peak non-gravitational acceleration      |
| **Outgassing Jet Offset**         | $8^\\circ$                                                                     | Angle of jet misalignment                |
| **Cumulative Outgassing Drift**   | $236.36\\text{ m}$                                                             | Spatial deviation from pure gravity path |

---

## 🧮 Mathematical Formulation

### 1\. Perturbed Equations of Motion

The state vector $\\mathbf{y} = [x, y, z, v\_x, v\_y, v\_z]$ is integrated in heliocentric coordinates:

$$\\frac{d^2\\mathbf{r}}{dt^2} = -\\frac{\\mu\_\\odot}{r^3}\\mathbf{r} + \\mathbf{a}*{\\text{rad}} + \\mathbf{a}*{\\text{thrust}}$$

### 2\. Radiation &amp; Outgassing Acceleration

* Solar Radiation Acceleration ($\\mathbf{a}*{\\text{rad}}$): $$P*{\\text{rad}} = \\frac{S\_0}{c \\cdot r\_{\\text{AU}}^2} \\cdot A, \\quad \\mathbf{a}*{\\text{rad}} = \\frac{P*{\\text{rad}}}{M\_0} \\hat{\\mathbf{r}}$$
* Jet Outgassing Thrust ($\\mathbf{a}*{\\text{thrust}}$ for $r*{\\text{AU}} &lt; 2.5$): $$\\mathbf{a}*{\\text{thrust}} = 0.2 \\cdot a*{\\text{rad}} \\left[ \\cos(8^\\circ)\\hat{\\mathbf{r}} + \\sin(8^\\circ)\\hat{\\mathbf{t}} \\right]$$

---

## 💻 Repository Contents

* `3iatlas_simulation_code-v2.py`: Main execution script containing the `RK45` solver, physical constants, differential equations, and visualization pipeline.
* `README.md`: Project documentation and mathematical guide.

### Auto-Generated Figures (Saved upon execution):

1. `trajectory_3iatlas.png`: 2D Heliocentric trajectory plotted against Earth ($1.0\\text{ AU}$) and Mars ($1.52\\text{ AU}$) orbits.
2. `velocity_3iatlas.png`: Velocity magnitude evolution vs. time compared to $v\_\\infty = 64.04\\text{ km/s}$.
3. `energy_3iatlas.png`: Specific orbital energy $\\epsilon = \\frac{1}{2}v^2 - \\frac{\\mu\_\\odot}{r}$ over the 138-day integration span.

---

## 🚀 Getting Started &amp; Execution

### Prerequisites

Make sure you have Python 3.8+ installed along with the required scientific packages:

```
pip install numpy scipy matplotlib

```

### Running the Code

Run the simulation script directly from your terminal or Google Colab:

```
python 3iatlas_simulation_code-v2.py

```

---

## ✍️ Author &amp; Citation

**Rene Sagal Andrade**  
*ORCID*: [0009-0003-5574-5348](https://www.google.com/url?sa=E&amp;q=https%3A%2F%2Forcid.org%2F0009-0003-5574-5348)  
*Email*: sagal233rene@outlook.com

If you use or reference this simulation code or the 3iatlas platform in your research, please cite the accompanying paper:

&gt; Sagal Andrade, R. (2026). *High-Precision Numerical Validation of Hyperbolic Trajectories: Assessing Non-Gravitational Forces in Interstellar Comet 3I/ATLAS*. Submitted to *Icarus* (Elsevier).
