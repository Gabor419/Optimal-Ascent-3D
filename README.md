# 3D 3-Optimal Ascent Trajectory Optimization (ICLOCS2 + IPOPT)

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: Proprietary](https://img.shields.io/badge/License-Proprietary-red.svg)](LICENSE)
[![Optimization: IPOPT](https://img.shields.io/badge/Solver-IPOPT%20%2F%20ICLOCS2-orange.svg)](https://imperialcollegelondon.github.io/ICLOCS/)

## Overview

This repository contains a MATLAB-based framework for solving the **3D 3-Degree-of-Freedom (3-DOF) Optimal Ascent Trajectory** problem for a multi-stage launch vehicle. 

This work represents an attempt to adapt and expand upon the problem described in the report for the *Launchers and Re-entry* course at Alma Mater Studiorum - UniBo, taught by Prof. Marco Sagliano. While the original report addressed the problem by solving a Two-Point Boundary Value Problem (TPBVP) using **Primer Vector Theory** (indirect method) via MATLAB's built-in functions, this script employs a **direct collocation method** via **ICLOCS2** and **IPOPT** to automatically optimize the trajectory subject to complex dynamic constraints and boundary conditions.

---

## Key Features

- **Direct Collocation Architecture:** Formulates the trajectory optimization as a Non-Linear Programming (NLP) problem, solved efficiently using IPOPT.
- **3D Rotating Earth Dynamics:** Full 3D equations of motion including gravity, Coriolis and centrifugal accelerations, aerodynamic drag, and spherical Earth kinematics.
- **Multi-Stage Modeling:** Explicit handling of 3 stages with continuous mass/thrust transitions utilizing smooth step functions ($\tanh$).
- **Verification Pipeline:** Includes built-in `ode45` simulation verification to validate the optimized control sequence against open-loop propagation.
- **Advanced Post-Processing:** Automated extraction of Cartesian costates (primer vector components), Hamiltonian evaluation, local (NEU), and global (ECEF) trajectory visualizations.

---

## Repository Structure

- `main.m` — Main script configuring vehicle data, executing the optimization solver, verifying via ODE simulation, and generating comprehensive analysis plots.
- `ascent.m` — Problem definition file containing boundaries, states, inputs, cost function definitions, and staging configurations for ICLOCS2.
- `ascent_dynamics.m` — 3-DOF equations of motion of the launch vehicle in spherical rotating coordinates.
- `settings_ascent.m` — Configuration file for ICLOCS2, containing transcription choices, meshing parameters, and customized IPOPT solver settings.

---

## Prerequisites & Dependencies

To run this code, you need:
1. **MATLAB** (tested on recent releases).
2. **ICLOCS2 Toolbox** (Imperial College Optimal Control Software v2). Ensure the `src` and `usr` directories are properly added to your MATLAB path.
3. **IPOPT** NLP solver (integrated via ICLOCS2).

---

## Usage

1. Open MATLAB and clone/add this repository to your workspace.
2. Update the paths in `main.m` to point to your local installation of `ICLOCS2`.
3. Run the `main.m` script:
   ```matlab
   run('main.m');
   ```

---

## Note on Geographic & Launch Site Selection

The choice of the **Semnan Space Centre** as the launch facility is strictly driven by physical, geometric, and dynamic problem requirements (such as latitude, inclination, and trajectory optimization mechanics). This code is developed solely for academic research and educational purposes, carrying no political or strategic affiliation.

---

## Author & License

* **Author:** Gabriele Borghini (2026)
* **Copyright:** © 2026 Gabriele Borghini. All Rights Reserved.
* **Proprietary & Confidential:** The source code contained herein is the property of the author. It may not be copied, modified, distributed, or disclosed to any third party without express written consent.