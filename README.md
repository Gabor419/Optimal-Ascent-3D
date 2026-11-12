# 3D 3S-Optimal Ascent Trajectory Optimization (ICLOCS2 + IPOPT)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#author--license)
[![Optimization: IPOPT](https://img.shields.io/badge/Solver-IPOPT%20%2F%20ICLOCS2-orange.svg)](https://imperialcollegelondon.github.io/ICLOCS/)

## Overview

This repository contains a MATLAB-based framework for solving the **3D 3-Degree-of-Freedom (3-DOF) Optimal Ascent Trajectory** problem for a three-stage (3S) launch vehicle.

This work represents an evolution and adaptation of a project originally developed in collaboration with Lisa Puttini for the *Launchers and Re-entry* course at **Alma Mater Studiorum - Università di Bologna**, taught by Prof. Marco Sagliano. While the original report addressed the problem by solving a Two-Point Boundary Value Problem (TPBVP) using **Primer Vector Theory** (indirect method) via MATLAB's built-in functions, this framework employs a **direct collocation method** with automatic mesh refinement via **ICLOCS2** and **IPOPT** to optimize the ascent trajectory subject to nonlinear path and terminal orbital insertion constraints.

---

## Key Features

- **Direct Collocation & Mesh Refinement:** Formulates the multi-stage trajectory optimization as a Non-Linear Programming (NLP) problem using Hermite-Simpson transcription and aggressive automatic mesh refinement (up to 10 iterations), solved via IPOPT with exact Hessian evaluation.
- **3D Rotating Earth Dynamics:** Full 3-DOF equations of motion over a rotating spherical Earth, accounting for inverse-square gravity, Coriolis and centrifugal accelerations, and exponential atmospheric drag.
- **Continuous Multi-Stage Modeling:** Explicit handling of 3 stages via a continuous mass state ($m_{\text{state}}$, varying solely due to propellant consumption) paired with smooth hyperbolic tangent ($\tanh$) step functions ($w = 1.5\text{ s}$) for thrust, mass flow rate, and structural mass-jettison transitions.
- **Right-Handed UEN Control Formulation:** Employs normalized local Cartesian thrust direction components ($\mathbf{u} = [u_U, u_E, u_N]$ in the right-handed **Up-East-North** frame) subject to unit-magnitude path bounds ($u_U^2 + u_E^2 + u_N^2 = 1$) and rate constraints ($\dot{\mathbf{u}} \in [-0.05, 0.05]\text{ s}^{-1}$), avoiding trigonometric steering singularities at vertical liftoff.
- **Exact Terminal Orbital Insertion:** Minimizes final burn time $t_f$ (equivalent to maximizing final injected mass during constant-thrust Stage 3 flight) while enforcing target altitude ($h_f = 300\text{ km}$), zero flight-path angle ($\gamma_f = 0^\circ$), target orbital inclination ($i = 35.66^\circ$), and circular inertial orbital energy.
- **Rigorous Verification Pipeline:** Includes built-in open-loop propagation via MATLAB's `ode45` integrator with strict tolerances ($\text{RelTol} = \text{AbsTol} = 10^{-10}$) and unit-norm control re-projection to validate the optimized trajectory against continuous dynamics.
- **Advanced Post-Processing:** Automated mapping of spherical Lagrange multipliers to ECEF Cartesian position and velocity costates (Primer Vector analysis via the right-handed UEN rotation matrix with $\det(R) = +1$), Hamiltonian evaluation, local right-handed curvilinear (**ENU**: East-North-Up) trajectory plots, and 3D global (**ECEF**) visualizations comparing the achieved vs. theoretical target orbit.

---

## Repository Structure

- `main.m` — Main orchestration script configuring vehicle parameters, executing the ICLOCS2/IPOPT solver with mesh refinement, performing open-loop `ode45` verification, and generating kinematic, costate, local ENU, and global 3D ECEF orbit insertion plots.
- `ascent.m` — Optimal Control Problem (OCP) definition file specifying state/control/rate bounds, piecewise-linear initial guess generation, staging timelines, minimum-time/maximum-mass objective, and terminal orbital boundary conditions (`b_unscaled`).
- `ascent_dynamics.m` — Vectorized 3-DOF equations of motion in spherical rotating coordinates driven by right-handed UEN local Cartesian thrust components and smooth $\tanh$ staging transitions.
- `settings_ascent.m` — Configuration file for ICLOCS2 (Hermite-Simpson discretization, automatic mesh refinement) and tuned IPOPT solver tolerances ($10^{-7}$).

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

The choice of the **Semnan Space Centre** ($\phi_0 = 35.23^\circ\text{ N}$, $\theta_0 = 53.92^\circ\text{ E}$, launch azimuth $\chi_0 = 84.08^\circ$) as the launch facility is strictly driven by physical, geometric, and dynamic problem requirements (such as launch latitude, target orbital inclination, and trajectory optimization mechanics). This code is developed solely for academic research and educational purposes, carrying no political or strategic affiliation.

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
- **License (MIT):** Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, and sublicense copies of the Software, subject to the inclusion of the original copyright and permission notice in all copies or substantial portions of the Software.
- **Disclaimer:** This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of accuracy, numerical stability, or fitness for a particular purpose.