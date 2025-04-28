% =========================================================
% Internal Ballistics Simulator - Refactored Version
% =========================================================
% core/runSimulation.m
% Function responsible for:
% 1. Setting ODE solver options (including the event function).
% 2. Defining initial conditions for state variables.
% 3. Preparing parameters to pass to the ODE function and event.
% 4. Calling the ODE solver (e.g., ode45).
% 5. Extracting raw results from the solver.
% 6. Calculating derived quantities post-simulation (pressure, gas density, etc.).
% 7. Recalculating torques over time from results (if necessary/desired).
% 8. Returning complete results in an organized struct.
% =========================================================

function results = runSimulation(parameters)
% INPUT:
%   parameters: Complete struct containing all necessary parameters
%               (output from loadCaseParameters).
% OUTPUT:
%   results: Struct containing simulation time vectors.

disp('Initializing ODE solver...');

% --- Validate required parameters ---
required_params_for_run = {'barrelLength_m', 'initialPropellantMass_m', 'initialTemperature_K', 'maxSafetyTime_s', ...
                           'projArea_Ab', 'initialFreeVolume_V0', 'covolume_b', 'specificGasConstant_R', ...
                           'twistRate_rad_m', 'engravingTorque_Nm', 'engravingEndPosition_m', ...
                           'rotationalFrictionCoeff_muRot', 'boreDiameter_Db', 'shotStartPressure_Pa', ...
                           'projMomentOfInertia_Ip'}; % Added more checks
missing_fields = required_params_for_run(~isfield(parameters, required_params_for_run));
if ~isempty(missing_fields)
    error('runSimulation: Missing required fields in parameters struct: %s', strjoin(missing_fields, ', '));
end

% --- ODE Solver Settings ---
options = odeset(...
    'RelTol', 1e-5, ...
    'AbsTol', 1e-7, ... % Adjust tolerances as needed
    'NonNegative', [1, 3, 4, 5, 6, 7], ... % Indices for non-negative states: m_prop_rem, x_p, v_p, omega, W_bore_res, Q_loss
    'Events', @(t,y) barrelExitEvent(t, y, parameters.barrelLength_m) ... % Use correct field name
);
fprintf('ODE options set. Event: barrelExitEvent @ x = %.4f m\n', parameters.barrelLength_m);

% --- Initial Conditions ---
% y = [m_prop_rem, T_gas, x_p, v_p, omega, W_bore_res, Q_loss]
initialConditions = [
    parameters.initialPropellantMass_m, ...   % Use correct field name
    parameters.initialTemperature_K, ...       % Use correct field name
    0, ...                                   % Initial projectile position [m]
    0, ...                                   % Initial projectile velocity [m/s]
    0, ...                                   % Initial angular velocity [rad/s]
    0, ...                                   % Initial accumulated bore resistance work [J]
    0                                        % Initial accumulated heat loss [J]
];
disp('Initial conditions defined:');
disp(initialConditions);

% --- Time Span ---
timeSpan = [0, parameters.maxSafetyTime_s]; % Use correct field name [s]
fprintf('Time interval: [%.1f, %.4f] s\n', timeSpan(1), timeSpan(2));

% --- Call ODE Solver ---
% The ODE function 'simulationOdes' receives the 'parameters' struct directly.
t = []; solution = []; eventTime = []; eventState = []; eventIndex = []; % Initialize ODE outputs

try
    disp('Starting ODE integration with ode45...');
    % Ensure 'simulationOdes' is the correct name of your ODE function
    [t, solution, eventTime, eventState, eventIndex] = ode45(@(t,y) simulationOdes(t, y, parameters), timeSpan, initialConditions, options);
    disp('ODE integration finished.');

    % Handle simulation end (event or max time)
    if ~isempty(eventIndex)
        fprintf('Barrel exit event detected at t = %.6f s (Index %d)\n', eventTime(end), eventIndex(end));
        % Optional: Truncate results at event time
        % t = t(1:eventIndex(end));
        % solution = solution(1:eventIndex(end),:);
        % disp('Results truncated to event time.');
        % Note: Currently keeping all points returned by ode45
    elseif ~isempty(t)
         warning('Simulation stopped at maximum safety time t_max_safety (%.4f s) BEFORE reaching barrel length (%.4f m).', t(end), parameters.barrelLength_m); % Use correct field name
    else
         error('ODE solver ode45 did not produce results.');
    end

catch ME_ode
    fprintf(2, 'Error during ode45 call or simulationOdes function: %s\n', ME_ode.message);
     if ~isempty(ME_ode.stack)
        fprintf(2, 'File: %s, Line: %d\n', ME_ode.stack(1).file, ME_ode.stack(1).line);
     end
    results = []; % Return empty on ODE failure
    return;
end

% --- Extract ODE Results and Post-Simulation Calculations ---
% Ensure the solution has the expected number of state variables from simulationOdes
numStateVarsExpected = 7; % [m_prop_rem; T_gas; x_p; v_p; omega; W_bore_res; Q_loss]
if isempty(t) || size(solution, 2) ~= numStateVarsExpected
     error('Output from ode45 is invalid or has incorrect number of columns (expected %d).', numStateVarsExpected);
end

disp('Extracting ODE results and calculating derived quantities...');
results = struct();
% --- Store results using consistent names expected by downstream functions ---
results.timeS                     = t(:); % Ensure column vector
results.remainingPropellantMassKg = solution(:, 1);
results.gasTemperatureK           = solution(:, 2);
results.projectilePositionM       = solution(:, 3);
results.projectileVelocityMps     = solution(:, 4);
results.angularVelocityRadps      = solution(:, 5);
results.frictionWorkJ             = solution(:, 6); % NOTE: This state variable is named W_bore_res in the ODE, represents work against bore resistance
results.heatLossJ                 = solution(:, 7);

% Calculate derived quantities (Pressure, Volume, Gas Density, Gas Mass)
results.gasMassKg = parameters.initialPropellantMass_m - results.remainingPropellantMassKg; % Use correct field name
results.gasMassKg = max(results.gasMassKg, 1e-12); % Avoid division by zero or log of zero

results.instantaneousVolumeM3 = parameters.initialFreeVolume_V0 + parameters.projArea_Ab .* results.projectilePositionM; % Use correct field names
results.instantaneousVolumeM3 = max(results.instantaneousVolumeM3, 1e-12);

effectiveVolume = results.instantaneousVolumeM3 - results.gasMassKg .* parameters.covolume_b; % Use correct field name
effectiveVolume = max(effectiveVolume, 1e-12); % Effective volume cannot be <= 0

results.gasPressurePa = (results.gasMassKg .* parameters.specificGasConstant_R .* results.gasTemperatureK) ./ effectiveVolume; % Use correct field names [Pa]
results.gasPressurePa = max(results.gasPressurePa, 1e-3); % Avoid negative or too low pressures

results.gasDensityKgm3 = results.gasMassKg ./ results.instantaneousVolumeM3; % Average gas density [kg/m^3]
results.gasDensityKgm3 = max(results.gasDensityKgm3, 1e-6);

disp('Derived quantities (gasPressure, gasDensity, instantaneousVolume, gasMass) calculated.');

% --- Recalculate Torques ---
% This block recalculates torques using simulation results.
% Useful for having exact values at each saved time step.
fprintf('Recalculating torques from ODE solution...\n');
numSteps = length(results.timeS);
results.riflingTorqueNm    = zeros(numSteps, 1); % Rifling torque [N*m]
results.engravingTorqueNm  = zeros(numSteps, 1); % Engraving torque [N*m]
results.frictionTorqueNm   = zeros(numSteps, 1); % Rotational friction torque [N*m]
results.netTorqueNm        = zeros(numSteps, 1); % Net torque [N*m]

% Parameters needed for torque calculation (use names from input 'parameters' struct)
projArea_Ab         = parameters.projArea_Ab;           % Use correct field name
twistRate_rad_m     = parameters.twistRate_rad_m;       % Use correct field name
engravingTorque_Nm  = parameters.engravingTorque_Nm;    % Use correct field name
engravingEndPos_m   = parameters.engravingEndPosition_m;% Use correct field name
rotFrictionCoeff    = parameters.rotationalFrictionCoeff_muRot; % Use correct field name
boreDiameter_Db     = parameters.boreDiameter_Db;       % Use correct field name
shotStartPressurePa = parameters.shotStartPressure_Pa;  % Use correct field name
% inertiaKgM2         = parameters.projMomentOfInertia_Ip; % Not directly used here, but available

boreRadiusM = boreDiameter_Db / 2.0; % Use correct field name
% Precalculate constant factors if possible
% Simplified form used in original ODEs: Fp * tan(alpha) * r = Fp * (twist * r) * r = Fp * twist * r^2
riflingFactor = twistRate_rad_m * boreRadiusM^2; % Use correct field name
rotFricFactor = rotFrictionCoeff * boreRadiusM;   % Use correct field name

for i = 1:numSteps
    currentPressurePa = results.gasPressurePa(i);
    currentPositionM  = results.projectilePositionM(i);

    % Apply activation logic (e.g., Pressure > Shot Start Pressure)
    if currentPressurePa > shotStartPressurePa && currentPositionM >= 0 % Ensure pressure and potential motion

        forcePressureI = currentPressurePa * projArea_Ab; % Use correct field name

        % Rifling Torque
        results.riflingTorqueNm(i) = forcePressureI * riflingFactor;

        % Engraving Torque
        if currentPositionM < engravingEndPos_m % Use correct field name
            results.engravingTorqueNm(i) = engravingTorque_Nm; % Use correct field name
        else
            results.engravingTorqueNm(i) = 0;
        end

        % Rotational Friction Torque (Assuming Normal Force proportional to Pressure Force)
        results.frictionTorqueNm(i) = forcePressureI * rotFricFactor; % Use pressure force as proxy

        % Net Torque (Accelerating - Resisting)
        torqueResistI = max(0, results.engravingTorqueNm(i) + results.frictionTorqueNm(i)); % Ensure non-negative
        results.netTorqueNm(i) = results.riflingTorqueNm(i) - torqueResistI;
        % Note: This net torque should ideally match I_p * d(omega)/dt calculated by the ODEs.
        % Slight differences can exist due to numerical integration.

    else
        % If pressure too low or projectile potentially moving backward (not typical here), torques are zero
        results.riflingTorqueNm(i)   = 0;
        results.engravingTorqueNm(i) = 0;
        results.frictionTorqueNm(i)  = 0;
        results.netTorqueNm(i)       = 0;
    end
end
disp('Torque recalculation completed.');

% --- Add event information to results struct ---
results.eventTimeS   = eventTime; % Event time(s)
results.eventState   = eventState;% State(s) at event(s)
results.eventIndex = eventIndex;% Event index/indices

disp('Final results struct assembled.');

end % End of function runSimulation