% =========================================================
% Internal Ballistics Simulator - Refactored Version
% =========================================================
% postprocessing/displayInputParameters.m
% Function dedicated to the formatted display of simulation
% input parameters on the console.
% =========================================================

function displayInputParameters(simulationParameters)
% INPUT:
%   simulationParameters: Complete struct containing input parameters
%                         (output of loadCaseParameters.m).

disp('--- Simulation Input Parameters ---');

% --- Case and File Identification Block ---
fprintf('IDENTIFICATION:\n');
fprintf('  Loaded Case Name          : %s\n', simulationParameters.loadedCaseName); % Assumes refactored name
if isfield(simulationParameters, 'benchmarkName') % Assumes refactored name
    fprintf('  Benchmark Reference       : %s\n', simulationParameters.benchmarkName);
end
if isfield(simulationParameters, 'propellantName') % Assumes refactored name
    fprintf('  Propellant Name/Type      : %s\n', simulationParameters.propellantName);
end
fprintf('  Case Config File          : %s\n', simulationParameters.caseConfigFile); % Assumes refactored name
fprintf('  Powder Data File          : %s\n', simulationParameters.powderDataFile); % Use powderDataFile
fprintf('\n');

% --- Geometric and Mass Properties Block (SI Units) ---
fprintf('GEOMETRIC AND MASS PROPERTIES (SI Units):\n');
fprintf('  Projectile Mass (m_p)     : %.6f kg (%.1f gr)\n', simulationParameters.projectileMassKg, simulationParameters.projectileWeightGr);
fprintf('  Projectile Base Area (A_b): %.3e m^2\n', simulationParameters.baseAreaM2);
fprintf('  Bore Diameter (D_bore)    : %.4f mm (%.3f inches)\n', simulationParameters.boreDiameterM * 1000, simulationParameters.projectileDiameterIn);
fprintf('  Axial Inertia (I_p)       : %.4e kg*m^2\n', simulationParameters.projectileInertiaKgM2);
fprintf('  Barrel Length (L_b)       : %.4f m\n', simulationParameters.barrelLength);
fprintf('  Propellant Mass (m_i)     : %.6f kg (%.1f gr)\n', simulationParameters.initialPropellantMassKg, simulationParameters.propellantMassGr);
fprintf('  Initial Free Vol. (V0)    : %.3e m^3 (%.2f cm^3)\n', simulationParameters.initialFreeVolumeM3, simulationParameters.initialFreeVolumeM3 * 1e6);
fprintf('  Propellant Solid Dens.(rho_i): %.1f kg/m^3\n', simulationParameters.propellantDensityKgm3);
fprintf('  Initial Surf. Area (S0)   : %.3e m^2 (from sigma=%.3f m^2/kg)\n', simulationParameters.initialSurfaceAreaM2, simulationParameters.specificSurfaceAreaM2pkg); % Assuming refactored names
fprintf('\n');

% --- Thermodynamics and Gas Properties Block ---
fprintf('THERMODYNAMIC AND GAS PROPERTIES:\n');
fprintf('  Co-volume (b_i)           : %.4f m^3/kg\n', simulationParameters.coVolumeM3pkg);
fprintf('  Specific Energy (F_i)     : %.3e J/kg (Calculated/Provided)\n', simulationParameters.specificEnergyJpkg); % Assuming refactored name
fprintf('  Adiab. Flame Temp (T_0i)  : %.0f K\n', simulationParameters.adiabaticFlameTempK); % Assuming refactored name
fprintf('  Specific Heat Ratio (gamma): %.3f\n', simulationParameters.gammaGas); % Assuming refactored name
fprintf('  Specific Gas Const (R_gas): %.2f J/(kg*K) (Calculated)\n', simulationParameters.specificGasConstantJpkgK);
if isfield(simulationParameters, 'meanGasMolarMassKgmol') % Assuming refactored name
    fprintf('  Mean Gas Molar Mass(M_gas): %.4f kg/mol\n', simulationParameters.meanGasMolarMassKgmol);
end
fprintf('\n');

% --- Combustion Law and Form Function Block ---
fprintf('COMBUSTION LAW AND FORM FUNCTION:\n');
fprintf('  Combustion Rate Law (a * P^beta):\n');
fprintf('    Coefficient (a_i)       : %.3e [m/s * Pa^-(%.4f)]\n', simulationParameters.combustionCoeffA, simulationParameters.combustionExponentBeta); % Assuming refactored names
fprintf('    Exponent (beta_i)       : %.4f\n', simulationParameters.combustionExponentBeta);
fprintf('  Form Function (S/S0 = (1-f)(1+theta*f)):\n');
fprintf('    Form Parameter (theta)  : %.3f\n', simulationParameters.formFactorTheta); % Assuming refactored name
fprintf('\n');

% --- Motion, Friction, Engraving, Heat Loss Block ---
fprintf('MOTION, FRICTION, ENGRAVING, HEAT LOSS:\n');
fprintf('  Linear Motion Parameters:\n');
fprintf('    Bore Resistance (br)    : Defined by table (br_travel_m, br_pressure_Pa)\n'); % Indicate use of table
fprintf('    Shot Start Press. (P_ss): %.2f MPa\n', simulationParameters.shotStartPressurePa / 1e6);
fprintf('  Rotational Motion Parameters:\n');
fprintf('    Twist Rate (calculated) : %.4f rad/m (%.1f in/turn)\n', simulationParameters.twistRadPerMeter, simulationParameters.twistInPerTurn);
fprintf('    Engraving Torque (T_eng): %.2f N*m (ESTIMATE)\n', simulationParameters.engravingTorqueValue);
fprintf('    Engraving End (x_eng)   : %.4f m (ESTIMATE)\n', simulationParameters.engravingEndPositionM);
fprintf('    Rot. Fric. Coeff.(mu_r) : %.4f (ESTIMATE)\n', simulationParameters.rotationalFrictionCoefficient);
fprintf('  Heat Loss Parameters:\n');
fprintf('    Convective Coeff. (h)   : %.1f W/(m^2*K) (ESTIMATE)\n', simulationParameters.convectionCoefficient);
fprintf('\n');

% --- Initial Conditions and Simulation Control Block ---
fprintf('INITIAL CONDITIONS AND SIMULATION CONTROL:\n');
fprintf('  Initial Sys. Temp (T0)    : %.2f K\n', simulationParameters.initialTemperatureK);
fprintf('  Max Simulation Time(t_max): %.4f s\n', simulationParameters.maxSafetyTimeS);
fprintf('\n');

disp('--------------------------------------------');

end % End of function displayInputParameters