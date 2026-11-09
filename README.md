# 3D 3S-Optimal Ascent Trajectory Optimization (ICLOCS2 + IPOPT)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: Proprietary](https://img.shields.io/badge/License-Proprietary-red.svg)](#author--license)
[![Optimization: IPOPT](https://img.shields.io/badge/Solver-IPOPT%20%2F%20ICLOCS2-orange.svg)](https://imperialcollegelondon.github.io/ICLOCS/)

## Overview

This repository contains a MATLAB-based framework for solving the **3D 3-Degree-of-Freedom (3-DOF) Optimal Ascent Trajectory** problem for a three-stage (3S) launch vehicle.

This work represents an evolution and adaptation of a project originally developed in collaboration with Lisa Puttini for the *Launchers and Re-entry* course at **Alma Mater Studiorum - Università di Bologna**, taught by Prof. Marco Sagliano. While the original report addressed the problem by solving a Two-Point Boundary Value Problem (TPBVP) using **Primer Vector Theory** (indirect method) via MATLAB's built-in functions, this framework employs a **direct collocation method** with automatic mesh refinement via **ICLOCS2** and **IPOPT** to optimize the ascent trajectory subject to nonlinear path and terminal orbital insertion constraints.

---

## Key Features

- **Direct Collocation & Mesh Refinement:** Formulates the multi-stage trajectory optimization as a Non-Linear Programming (NLP) problem using Hermite-Simpson transcription and automatic mesh refinement, solved via IPOPT with exact Hessian evaluation.
- **3D Rotating Earth Dynamics:** Full 3-DOF equations of motion over a rotating spherical Earth, accounting for inverse-square gravity, Coriolis and centrifugal accelerations, and exponential atmospheric drag.
- **Continuous Multi-Stage Modeling:** Explicit handling of 3 stages with continuous thrust, mass flow rate, and structural mass-jettison transitions modeled via smooth hyperbolic tangent ($\tanh$) step functions.
- **Singularity-Free Control Formulation:** Employs normalized local Cartesian thrust direction components ($\mathbf{u} = [u_x, u_y, u_z]$ in the Up-North-East frame) enforced by a unit-magnitude path constraint ($\|\mathbf{u}\|^2 = 1$), avoiding trigonometric steering singularities at vertical liftoff.
- **Rigorous Verification Pipeline:** Includes built-in open-loop propagation via MATLAB's `ode45` integrator with strict tolerances ($10^{-10}$) to validate the optimized control law against continuous dynamics.
- **Advanced Post-Processing:** Automated mapping of spherical Lagrange multipliers to ECEF Cartesian costates (Primer Vector analysis), Hamiltonian evaluation, local curvilinear (NEU) plots, and 3D global (ECEF) orbit insertion visualizations.

---

## Repository Structure

- `main.m` — Main orchestration script configuring vehicle parameters, executing the ICLOCS2/IPOPT solver, performing open-loop `ode45` verification, and generating kinematic, costate, and 3D trajectory plots.
- `ascent.m` — Optimal Control Problem (OCP) definition file specifying state/control bounds, initial guess generation, staging timelines, minimum-time/maximum-mass objective, and terminal orbital constraints (target altitude, inclination, and circular orbital energy).
- `ascent_dynamics.m` — Vectorized 3-DOF equations of motion in spherical rotating coordinates driven by local Cartesian thrust components and smooth staging transitions.
- `settings_ascent.m` — Configuration file for ICLOCS2 (Hermite-Simpson discretization, automatic mesh refinement) and bespoke IPOPT solver tolerances.

---

## Prerequisites & Dependencies

To run this code, you need:
1. **MATLAB** (R2021b or later recommended, as the plotting suite uses the `theme` command).
2. **ICLOCS2 Toolbox** (Imperial College London Optimal Control Software v2). Ensure the `src` and `usr` directories are properly added to your MATLAB path.
3. **IPOPT** NLP solver (interfaced via MEX files within ICLOCS2).
4. **ADiGator** *(optional)* — Only required if switching derivative generation from `'numeric'` to `'adigator'` in `settings_ascent.m`.

---

## Usage

1. Clone or download this repository and add it to your MATLAB workspace.
2. Update the path configurations at the top of `main.m` (and `settings_ascent.m` if using ADiGator) to point to your local installation of `ICLOCS2`.
3. Run the `main.m` script:
   ```matlab
   run('main.m');
   ```

---

## Note on Geographic & Launch Site Selection

The choice of the **Semnan Space Centre** as the launch facility is strictly driven by physical, geometric, and dynamic problem requirements (such as launch latitude, target orbital inclination, and trajectory optimization mechanics). This code is developed solely for academic research and educational purposes, carrying no political or strategic affiliation.

---

## Academic Context & Acknowledgments

This code is developed solely for academic research, advanced astrodynamics studies, and educational purposes at **Alma Mater Studiorum - Università di Bologna**.

This framework interfaces with third-party open-source optimization tools, which remain the property of their respective authors under their own licensing terms:
- **ICLOCS2:** Y. Nie, O. Faqir, and E. C. Kerrigan (2018), *"ICLOCS2: Try this Optimal Control Problem Solver Before you Try the Rest"*, Imperial College London (MIT License).
- **IPOPT:** A. Wächter and L. T. Biegler (2006), *"On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming"*, Mathematical Programming (EPL License).

---

## Author & License

- **Author:** Gabriele Borghini (2026)
- **Copyright:** © 2026 Gabriele Borghini. All Rights Reserved.
- **Usage Rights:** The source code, mathematical formulations, and custom scripts contained in this repository are the intellectual property of the author. They may not be copied, modified, redistributed, or used for commercial purposes without the express written consent of the author.
- **Disclaimer:** This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of accuracy, numerical stability, or fitness for a particular purpose.
