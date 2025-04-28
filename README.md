# TORSHOT: Internal Ballistics Simulator in MATLAB

**Development and Analysis of a Lumped Parameter Model for Predicting Small Caliber Ammunition Performance**

*Author: Francesco Torsello*
*Institution: Politecnico di Milano* [cite: 2]

---

## Overview

TORSHOT is a computational simulator for the internal ballistics of small-caliber ammunition, developed in the MATLAB environment[cite: 1, 5]. It aims to predict key performance metrics, such as chamber pressure profiles and projectile muzzle velocity, by simulating the complex thermo-fluid dynamic phenomena occurring inside a firearm during firing[cite: 3, 4, 5]. The simulator serves as a valuable tool for parametric analysis and preliminary design of weapon-ammunition systems, complementing experimental testing[cite: 4, 9, 671].

The simulator employs a lumped-parameter (0D) model, assuming spatial uniformity of thermodynamic properties (pressure, temperature) within the control volume (chamber behind the projectile)[cite: 6, 182, 183, 184, 672].

## Core Features

* **0D Lumped Parameter Model:** Efficiently simulates the internal ballistic cycle[cite: 6, 182].
* **Core Physics Models:**
    * Propellant Combustion: Vieille's Law ($r = \alpha P_g^{\beta}$)[cite: 6, 95, 674].
    * Gas Thermodynamics: Noble-Able Equation of State ($P(V - m_g b) = m_g R_g T$)[cite: 6, 143, 675].
    * Projectile Dynamics: Equations of motion for translation and rotation, including engraving resistance, friction, and rifling effects[cite: 6, 676].
    * Heat Transfer: Convective heat loss model to the barrel ($h_{conv} \propto P_g^n$)[cite: 6, 677].
* **MATLAB Implementation:** Developed entirely in MATLAB for accessibility and leveraging its numerical capabilities[cite: 5, 7].
* **Modular Architecture:** Code organized into functional folders (`config`, `core`, `gui`, `postprocessing`, `utils`, `output`) for maintainability and scalability[cite: 251, 678].
* **Numerical Solver:** Utilizes MATLAB's `ode45` adaptive step Runge-Kutta solver with event detection for precise muzzle exit determination[cite: 7, 679].
* **Graphical User Interface (GUI):** An interactive GUI facilitates case setup (modular component selection: cartridge, bullet, powder, barrel), parameter input, simulation execution, and results visualization[cite: 7, 257, 680].
* **Post-Processing:** Advanced analysis tools including:
    * Detailed energy balance calculation[cite: 8, 681].
    * Barrel temperature increase estimation[cite: 8, 681].
    * Barrel stress analysis (based on Lamé equations and Tresca criterion)[cite: 468, 701, 702].
* **Parameter Optimization:** Includes a grid search module for calibrating uncertain model parameters (e.g., Vieille coefficients $\alpha, \beta$) against reference data[cite: 8, 260, 682].

## Model Overview

The simulation solves a system of coupled Ordinary Differential Equations (ODEs) describing the temporal evolution of the system state[cite: 6, 193]. Key state variables include remaining propellant mass, gas temperature, projectile position, linear velocity, and angular velocity[cite: 196, 197, 198, 199]. The core equations integrate Vieille's law for burn rate, the Noble-Able equation for gas state, Newton's second law for projectile motion (linear and angular), and simplified models for heat loss and bore resistance[cite: 674, 675, 676, 677].

## Software Structure

The codebase is organized as follows[cite: 251, 678]:

* `config/`: Contains configuration files for components (cartridges, bullets, powders, barrels) and default physics parameters[cite: 253, 254].
* `core/`: Houses the main simulation engine, ODE function definitions (`simulationOdes.m`), and the solver execution logic (`runSimulation.m`)[cite: 255, 256].
* `gui/`: Contains all files related to the graphical user interface[cite: 257].
* `postprocessing/`: Includes scripts for analyzing results, calculating energy balance, stress analysis, and parameter optimization[cite: 258, 259, 260].
* `utils/`: General utility functions (e.g., path management, safe field access)[cite: 261, 262].
* `output/`: Default directory for saving simulation results and summary files[cite: 264, 265].

## Getting Started

### Prerequisites

* MATLAB environment (specify version if known, requires standard toolboxes including ODE solvers).

### Usage

1.  **Launch:** Run the main GUI script (likely `ballisticSimulatorGUI.m` located in the `gui` folder)[cite: 257].
2.  **Configure:**
    * Use the dropdown menus in the GUI to select the desired cartridge, bullet, powder, and barrel definition files located in the `config` subfolders[cite: 269, 282]. Load these components.
    * Input simulation-specific parameters like propellant charge weight, barrel length, twist rate, initial temperature, etc., in the editable fields[cite: 272, 302, 303, 304].
    * Apply settings and calculate initial volumes using the appropriate button[cite: 274].
3.  **Run Simulation:** Execute the simulation using the "Run Simulation" button[cite: 275].
4.  **Analyze Results:** View plotted results in the GUI, check the command window for summaries, and explore detailed energy balance or stress analysis windows[cite: 277, 413, 414, 415, 416]. Saved results can be found in the `output` folder[cite: 264, 429, 430].

## Dependencies

* MATLAB (Version RXXXXx or later recommended)
* *(Add any specific MATLAB Toolboxes if required, e.g., Optimization Toolbox for grid search)*

## Limitations

* **0D Model:** Cannot capture spatial variations (pressure waves, temperature gradients, erosive burning)[cite: 704].
* **Predictive Accuracy:** Shows systematic deviations compared to experimental data (especially muzzle velocity), potentially due to model simplifications[cite: 638, 706]. Best suited for comparative analysis and trend evaluation unless carefully calibrated.
* **Parameter Sensitivity:** Results are highly sensitive to input parameters, particularly propellant properties ($\alpha, \beta, b, f$) and auxiliary models (heat transfer, friction), which require accurate calibration[cite: 698, 699, 700, 707].
* **Auxiliary Models:** Simplified models for heat loss and bore resistance/friction significantly impact results[cite: 650, 710, 711, 712].
* **Ignition Phase:** The complex ignition process is currently neglected[cite: 714].

## Future Work

Potential areas for future development include[cite: 716]:

* **Model Improvements:** Refine heat transfer and friction models, explore quasi-dimensional (1D) or two-phase flow approaches, model the ignition phase[cite: 717, 719, 720, 721, 722].
* **Functionality Extensions:** Integrate intermediate/external ballistics, expand component databases, add models for phenomena like barrel erosion[cite: 723, 724, 726].
* **Usability:** Enhance GUI, optimize code, improve documentation[cite: 727, 728, 729].
* **Validation:** Conduct large-scale validation against diverse experimental data, perform sensitivity/uncertainty quantification studies[cite: 730, 731, 732, 733].

