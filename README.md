# TORSHOT: Internal Ballistics Simulator in MATLAB

**Development and Analysis of a Lumped Parameter Model for Predicting Small Caliber Ammunition Performance**

## Overview

TORSHOT is a computational simulator for the internal ballistics of small-caliber ammunition, developed in the MATLAB environment. It aims to predict key performance metrics, such as chamber pressure profiles and projectile muzzle velocity, by simulating the complex thermo-fluid dynamic phenomena occurring inside a firearm during firing. The simulator serves as a valuable tool for parametric analysis and preliminary design of weapon-ammunition systems, complementing experimental testing.

The simulator employs a lumped-parameter (0D) model, assuming spatial uniformity of thermodynamic properties (pressure, temperature) within the control volume (chamber behind the projectile).

## Core Features

* **0D Lumped Parameter Model:** Efficiently simulates the internal ballistic cycle.
* **Core Physics Models:**
    * Propellant Combustion: Vieille's Law ($r = \alpha P_g^{\beta}$).
    * Gas Thermodynamics: Noble-Able Equation of State ($P(V - m_g b) = m_g R_g T$).
    * Projectile Dynamics: Equations of motion for translation and rotation, including engraving resistance, friction, and rifling effects.
    * Heat Transfer: Convective heat loss model to the barrel ($h_{conv} \propto P_g^n$).
* **MATLAB Implementation:** Developed entirely in MATLAB for accessibility and leveraging its numerical capabilities.
* **Modular Architecture:** Code organized into functional folders (`config`, `core`, `gui`, `postprocessing`, `utils`, `output`) for maintainability and scalability.
* **Numerical Solver:** Utilizes MATLAB's `ode45` adaptive step Runge-Kutta solver with event detection for precise muzzle exit determination.
* **Graphical User Interface (GUI):** An interactive GUI facilitates case setup (modular component selection: cartridge, bullet, powder, barrel), parameter input, simulation execution, and results visualization.
* **Post-Processing:** Advanced analysis tools including:
    * Detailed energy balance calculation.
    * Barrel temperature increase estimation.
    * Barrel stress analysis (based on Lamé equations and Tresca criterion).
* **Parameter Optimization:** Includes a grid search module for calibrating uncertain model parameters (e.g., Vieille coefficients $\alpha, \beta$) against reference data.

## Model Overview

The simulation solves a system of coupled Ordinary Differential Equations (ODEs) describing the temporal evolution of the system state. Key state variables include remaining propellant mass, gas temperature, projectile position, linear velocity, and angular velocity. The core equations integrate Vieille's law for burn rate, the Noble-Able equation for gas state, Newton's second law for projectile motion (linear and angular), and simplified models for heat loss and bore resistance.

## Software Structure

The codebase is organized as follows:

* `config/`: Contains configuration files for components (cartridges, bullets, powders, barrels) and default physics parameters.
* `core/`: Houses the main simulation engine, ODE function definitions (`simulationOdes.m`), and the solver execution logic (`runSimulation.m`).
* `gui/`: Contains all files related to the graphical user interface.
* `postprocessing/`: Includes scripts for analyzing results, calculating energy balance, stress analysis, and parameter optimization.
* `utils/`: General utility functions (e.g., path management, safe field access).
* `output/`: Default directory for saving simulation results and summary files.

## Getting Started

### Usage

1.  **Launch:** Run the main GUI script (likely `ballisticSimulatorGUI.m` located in the `gui` folder).
2.  **Configure:**
    * Use the dropdown menus in the GUI to select the desired cartridge, bullet, powder, and barrel definition files located in the `config` subfolders. Load these components.
    * Input simulation-specific parameters like propellant charge weight, barrel length, twist rate, initial temperature, etc., in the editable fields.
    * Apply settings and calculate initial volumes using the appropriate button.
3.  **Run Simulation:** Execute the simulation using the "Run Simulation" button.
4.  **Analyze Results:** View plotted results in the GUI, check the command window for summaries, and explore detailed energy balance or stress analysis windows. Saved results can be found in the `output` folder.

## Limitations

* **0D Model:** Cannot capture spatial variations (pressure waves, temperature gradients, erosive burning).
* **Predictive Accuracy:** Shows systematic deviations compared to experimental data (especially muzzle velocity), potentially due to model simplifications. Best suited for comparative analysis and trend evaluation unless carefully calibrated.
* **Parameter Sensitivity:** Results are highly sensitive to input parameters, particularly propellant properties ($\alpha, \beta, b, f$) and auxiliary models (heat transfer, friction), which require accurate calibration.
* **Auxiliary Models:** Simplified models for heat loss and bore resistance/friction significantly impact results.
* **Ignition Phase:** The complex ignition process is currently neglected.

## Future Work

Potential areas for future development include:

* **Model Improvements:** Refine heat transfer and friction models, explore quasi-dimensional (1D) or two-phase flow approaches, model the ignition phase.
* **Functionality Extensions:** Integrate intermediate/external ballistics, expand component databases, add models for phenomena like barrel erosion.
* **Usability:** Enhance GUI, optimize code, improve documentation.
* **Validation:** Conduct large-scale validation against diverse experimental data, perform sensitivity/uncertainty quantification studies.
