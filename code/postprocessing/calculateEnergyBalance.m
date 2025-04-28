% =========================================================
% Internal Ballistics Simulator - Refactored Version
% =========================================================
% postprocessing/calculateEnergyBalance.m
% Performs energy balance, estimates barrel temperature rise, and calculates BC.
% Includes calculation of engraving work and Ballistic Coefficient (BC).
% CORRECTED: Uses parameter names consistent with loadCaseParameters output.
% =========================================================

function [balance, barrelTempIncreaseK] = calculateEnergyBalance(results, parameters, barrelParameters, formFactor)
% INPUTS:
%   results: Struct with simulation results (from runSimulation).
%            Required fields: timeS, projectileVelocityMps, angularVelocityRadps,
%                             remainingPropellantMassKg, heatLossJ, frictionWorkJ,
%                             frictionTorqueNm, projectilePositionM.
%   parameters: Struct with simulation input parameters (from loadCaseParameters).
%               Required fields: initialPropellantMass_m, impetus_F, projMass_m,
%                                projMomentOfInertia_Ip, engravingTorque_Nm, engravingEndPosition_m,
%                                projWeight_gr, projDiameter_in. % CORRECTED Names
%   barrelParameters: Struct with barrel parameters.
%                     Required fields: massKg, specificHeatJkgK.
%   formFactor: Projectile form factor (relative to G1 or other standard).
%               This value MUST be provided from external data.
%
% OUTPUTS:
%   balance: Struct containing calculated energies [J] and BC:
%            - EstimatedPowderEnergyJ, ProjectileLinearKeJ, ProjectileRotationalKeJ,
%            - ProjectileTotalKeJ, EstimatedGasKeJ, GasVelocityFactorUsed, TotalHeatLossJ, LinearFrictionWorkJ,
%            - RotationalFrictionWorkJ, TotalFrictionWorkJ, EngravingWorkJ, UnaccountedEnergyJ
%            - SectionalDensityLbPerIn2: Calculated sectional density [lb/in^2]
%            - BallisticCoefficient: Calculated BC (using formFactor)
%   barrelTempIncreaseK: Estimated average increase in barrel temperature [K]

% --- Input Validation (Basic) ---
% Validate results struct fields (using names matching runSimulation output)
results_required = {'timeS', 'projectileVelocityMps', 'angularVelocityRadps', ...
                    'remainingPropellantMassKg', 'heatLossJ', 'frictionWorkJ', ...
                    'frictionTorqueNm', 'projectilePositionM'};
missing_results = results_required(~isfield(results, results_required));
if ~isempty(missing_results)
    error('calculateEnergyBalance: The ''results'' structure is missing required fields: %s', strjoin(missing_results, ', '));
end

% Validate parameters struct fields (using names matching loadCaseParameters output)
parameters_required = {'initialPropellantMass_m', 'impetus_F', 'projMass_m', 'projMomentOfInertia_Ip', ...
                       'engravingTorque_Nm', 'engravingEndPosition_m', ...
                       'projWeight_gr', 'projDiameter_in'}; % CORRECTED Names
missing_params = parameters_required(~isfield(parameters, parameters_required));
if ~isempty(missing_params)
    error('calculateEnergyBalance: The ''parameters'' structure is missing required fields: %s', strjoin(missing_params, ', '));
end

if ~isstruct(barrelParameters) || ~isfield(barrelParameters, 'massKg') || ~isfield(barrelParameters, 'specificHeatJkgK')
    error('calculateEnergyBalance: The ''barrelParameters'' structure is missing ''massKg'' or ''specificHeatJkgK''.');
end
% Validate the formFactor input
% if nargin < 4 || ~isnumeric(formFactor) || isempty(formFactor) || formFactor <= 0
%     error('calculateEnergyBalance: Input ''formFactor'' is missing, invalid, or non-positive.');
% end
% Other validations
if barrelParameters.massKg <= 0 || barrelParameters.specificHeatJkgK <= 0
    error('calculateEnergyBalance: Barrel mass and specific heat must be positive.');
end
if isempty(results.timeS) || length(results.timeS) < 2
    error('calculateEnergyBalance: Time vector in results is invalid or too short.');
end

% --- Calculate Energy Components ---
balance = struct(); % Initialize output structure

% 1. Estimated Powder Energy (based on impetus and BURNED mass)
% Use correct parameter name: initialPropellantMass_m
burnedPropellantMassKg = parameters.initialPropellantMass_m - results.remainingPropellantMassKg(end);
burnedPropellantMassKg = max(0, burnedPropellantMassKg); % Ensure non-negative
% Use correct parameter name: impetus_F
balance.EstimatedPowderEnergyJ = burnedPropellantMassKg * parameters.impetus_F; % Use burned mass


% 2. Final Projectile Kinetic Energy
exitVelocity = results.projectileVelocityMps(end);
exitOmega    = results.angularVelocityRadps(end);
% Use correct parameter names: projMass_m, projMomentOfInertia_Ip
balance.ProjectileLinearKeJ     = 0.5 * parameters.projMass_m * exitVelocity^2;
balance.ProjectileRotationalKeJ = 0.5 * parameters.projMomentOfInertia_Ip * exitOmega^2;
balance.ProjectileTotalKeJ      = balance.ProjectileLinearKeJ + balance.ProjectileRotationalKeJ;

% 3. Estimated Muzzle Gas Kinetic Energy (Approximate)
gasVelocityFactor = 1.7; % Multiplicative factor gas velocity vs projectile velocity (empirical)
exitGasMass = burnedPropellantMassKg; % Mass of gas is the mass of burned propellant
estimatedExitGasVelocity = gasVelocityFactor * exitVelocity;
balance.EstimatedGasKeJ    = 0.5 * exitGasMass * estimatedExitGasVelocity^2;
balance.GasVelocityFactorUsed = gasVelocityFactor;

% 4. Total Heat Loss (Q_loss)
balance.TotalHeatLossJ = results.heatLossJ(end); % Access from results is correct

% 5. Total Friction Work (Linear + Rotational)
balance.LinearFrictionWorkJ = results.frictionWorkJ(end); % Access from results is correct
% Calculate Rotational Friction Work
timeVector              = results.timeS(:);
rotFrictionTorqueVector = results.frictionTorqueNm(:); % Access from results is correct
omegaVector             = results.angularVelocityRadps(:); % Access from results is correct
if length(timeVector) ~= length(rotFrictionTorqueVector) || length(timeVector) ~= length(omegaVector)
    error('calculateEnergyBalance: Vector dimensions for t, TorqueF, omega do not match.');
end
rotFrictionPower = rotFrictionTorqueVector .* omegaVector;
rotFrictionPower(isnan(rotFrictionPower) | isinf(rotFrictionPower)) = 0; % Handle potential NaN/Inf
balance.RotationalFrictionWorkJ = trapz(timeVector, rotFrictionPower);
balance.TotalFrictionWorkJ      = balance.LinearFrictionWorkJ + balance.RotationalFrictionWorkJ;

% 6. Engraving Work (W_engrave)
% Use correct parameter name: engravingEndPosition_m
engravingIndices = results.projectilePositionM < parameters.engravingEndPosition_m;
if any(engravingIndices)
    engravingTime   = timeVector(engravingIndices);
    engravingOmega  = omegaVector(engravingIndices);
    % Use correct parameter name: engravingTorque_Nm
    engravingPower  = parameters.engravingTorque_Nm * engravingOmega;
    engravingPower(isnan(engravingPower) | isinf(engravingPower)) = 0; % Handle potential NaN/Inf
    if length(engravingTime) > 1
        balance.EngravingWorkJ = trapz(engravingTime, engravingPower);
    else
        balance.EngravingWorkJ = 0; % Not enough points to integrate
    end
else
    balance.EngravingWorkJ = 0; % Engraving phase not captured or non-existent
end

% 7. Unaccounted Energy
% Add check for residual internal energy if needed (depends on model complexity, not strictly required for this balance check)
% Example placeholder: U_gas_res = 0;
U_gas_res = 0; % Assuming we don't calculate residual gas internal energy in this version
balance.UnaccountedEnergyJ = balance.EstimatedPowderEnergyJ - balance.ProjectileTotalKeJ ...
                           - balance.EstimatedGasKeJ - balance.TotalHeatLossJ ...
                           - balance.TotalFrictionWorkJ - balance.EngravingWorkJ - U_gas_res;

% --- Ballistic Coefficient (BC) Calculation ---
try
    % Use correct parameter names: projWeight_gr, projDiameter_in
    weightGr   = parameters.projWeight_gr;
    diameterIn = parameters.projDiameter_in;

    if diameterIn <= 0 || weightGr <= 0
        error('calculateEnergyBalance: Invalid projectile weight or diameter for BC calculation.');
    end

    % Calculate SD in lb/in^2
    KG_PER_GRAIN = 1 / 15432.35835;
    KG_TO_LB     = 2.20462;
    massKg       = weightGr * KG_PER_GRAIN;
    massLb       = massKg * KG_TO_LB;
    sectionalDensityLbPerIn2 = massLb / (diameterIn^2);

    % Calculate BC using provided form factor
    ballisticCoefficient = sectionalDensityLbPerIn2 / formFactor;

    % Save to balance struct
    balance.SectionalDensityLbPerIn2 = sectionalDensityLbPerIn2;
    balance.BallisticCoefficient     = ballisticCoefficient;

catch ME_bc
    warning('calculateEnergyBalance:ErrorCalculatingBC', ...
            'Error during BC calculation: %s. BC not calculated.', ME_bc.message);
    balance.SectionalDensityLbPerIn2 = NaN;
    balance.BallisticCoefficient     = NaN;
end
% --- END BC CALCULATION ---

% --- Calculate Barrel Temperature Increase (Method: Qloss + Total Friction Work + Engraving Work) ---
% --- Calculate Barrel Temperature Increase (Method: Qloss + Total Friction Work + Engraving Work) ---
depositedBarrelEnergyJ = balance.TotalHeatLossJ + balance.TotalFrictionWorkJ + balance.EngravingWorkJ;

% --- DEBUGGING LOG STATUS ADDED ---
logStatus('>> DEBUG calculateEnergyBalance (Barrel Temp Calc):'); % Usa logStatus
logStatus(sprintf('   TotalHeatLossJ         = %.4f J', balance.TotalHeatLossJ)); % Usa logStatus
% Display Linear/Rotational breakdown if available, otherwise just TotalFrictionWorkJ
if isfield(balance, 'LinearFrictionWorkJ') && isfield(balance, 'RotationalFrictionWorkJ')
    logStatus(sprintf('   TotalFrictionWorkJ     = %.4f J (Linear: %.4f + Rotational: %.4f)', balance.TotalFrictionWorkJ, balance.LinearFrictionWorkJ, balance.RotationalFrictionWorkJ)); % Usa logStatus
else
    logStatus(sprintf('   TotalFrictionWorkJ     = %.4f J', balance.TotalFrictionWorkJ)); % Usa logStatus
end
logStatus(sprintf('   EngravingWorkJ         = %.4f J', balance.EngravingWorkJ)); % Usa logStatus
logStatus('   ----------------------------------'); % Usa logStatus
logStatus(sprintf('   Deposited Barrel Energy = %.4f J', depositedBarrelEnergyJ)); % Usa logStatus
logStatus(sprintf('   Barrel Mass             = %.4f kg', barrelParameters.massKg)); % Usa logStatus
logStatus(sprintf('   Barrel Specific Heat    = %.4f J/kgK', barrelParameters.specificHeatJkgK)); % Usa logStatus
logStatus(sprintf('   Denominator (Mass * Cp) = %.4f', (barrelParameters.massKg * barrelParameters.specificHeatJkgK))); % Usa logStatus
% --- END DEBUGGING LOG STATUS ---

barrelTempIncreaseK    = depositedBarrelEnergyJ / (barrelParameters.massKg * barrelParameters.specificHeatJkgK); % [K]

% Add another display right after calculation to confirm result
logStatus(sprintf('   Calculated barrelTempIncreaseK = %.4f K', barrelTempIncreaseK)); % Usa logStatus
logStatus('<< END DEBUG calculateEnergyBalance'); % Usa logStatus

end % End of function calculateEnergyBala