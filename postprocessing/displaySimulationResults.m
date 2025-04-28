% =========================================================
% Internal Ballistics Simulator - Refactored Version
% =========================================================
% postprocessing/displaySimulationResults.m
% Function dedicated to generating plots of the main
% simulation results.
% CORRECTED: Uses parameter names consistent with loadCaseParameters output.
% =========================================================

function displaySimulationResults(simulationResults, parameters)
% INPUTS:
%   simulationResults: Struct containing simulation results
%                      (output of runSimulation.m).
%   parameters: Struct containing simulation input parameters
%               (output of loadCaseParameters.m).

disp('Generating results plots...');

% --- Validate inputs ---
if ~isstruct(simulationResults) || ~isfield(simulationResults, 'timeS') || isempty(simulationResults.timeS)
    error('displaySimulationResults: Invalid or empty simulationResults struct provided.');
end
if ~isstruct(parameters) || ~isfield(parameters, 'loadedCaseName') % Basic check
    error('displaySimulationResults: Invalid or incomplete parameters struct provided.');
end
required_results_fields = {'gasPressurePa', 'projectilePositionM', 'projectileVelocityMps', ...
                           'angularVelocityRadps', 'remainingPropellantMassKg', 'gasMassKg', ...
                           'gasTemperatureK', 'frictionWorkJ', 'heatLossJ', 'riflingTorqueNm', ...
                           'engravingTorqueNm', 'frictionTorqueNm', 'netTorqueNm'};
missing_results = required_results_fields(~isfield(simulationResults, required_results_fields));
if ~isempty(missing_results)
    error('displaySimulationResults: Missing required fields in simulationResults: %s', strjoin(missing_results, ', '));
end
required_params_fields = {'loadedCaseName', 'barrelLength_m', 'convectiveHeatCoeff_h', ...
                          'boreDiameter_Db', 'initialTemperature_K'}; % Check required parameters
missing_params = required_params_fields(~isfield(parameters, required_params_fields));
if ~isempty(missing_params)
    error('displaySimulationResults: Missing required fields in parameters: %s', strjoin(missing_params, ', '));
end


% --- Extract Data for Plotting ---
timeS = simulationResults.timeS;
% Validate time data
if isempty(timeS) || length(timeS) <= 1
    warning('displaySimulationResults: No valid time data to plot.');
    return;
end

gasPressurePa            = simulationResults.gasPressurePa;
projectilePositionM      = simulationResults.projectilePositionM;
projectileVelocityMps    = simulationResults.projectileVelocityMps;
angularVelocityRadps     = simulationResults.angularVelocityRadps;
remainingPropellantMassKg= simulationResults.remainingPropellantMassKg;
gasMassKg                = simulationResults.gasMassKg;
gasTemperatureK          = simulationResults.gasTemperatureK;
frictionWorkJ            = simulationResults.frictionWorkJ; % Bore resistance work
heatLossJ                = simulationResults.heatLossJ;
riflingTorqueNm          = simulationResults.riflingTorqueNm;
engravingTorqueNm        = simulationResults.engravingTorqueNm;
frictionTorqueNm         = simulationResults.frictionTorqueNm;
netTorqueNm              = simulationResults.netTorqueNm;

% --- Preliminary Calculations for Plots ---
timeMs        = timeS * 1000;       % Time in milliseconds
pressureMPa   = gasPressurePa / 1e6; % Pressure in MPa
positionCm    = projectilePositionM * 100; % Position in cm
omegaRpm      = angularVelocityRadps * (60 / (2 * pi)); % Omega in RPM
propellantMassG= remainingPropellantMassKg * 1000; % Masses in grams
gasMassG      = gasMassKg * 1000;
frictionWorkKJ= frictionWorkJ / 1000; % Work/Heat in kJ (Bore resistance work)
heatLossKJ    = heatLossJ / 1000;

% Calculate Heat Loss Rate (dQ/dt) a posteriori (for visualization)
heatLossRateW = zeros(size(timeS));
% --- Use CORRECTED parameter names ---
convectionCoefficient = parameters.convectiveHeatCoeff_h;
boreDiameterM         = parameters.boreDiameter_Db;
initialBarrelTempK    = parameters.initialTemperature_K;
% --- End CORRECTED parameter names ---
estimatedChamberAreaM2= 1e-3; % USE THE SAME ESTIMATE AS IN THE ODEs! (Or pass/recalculate)

for i = 1:length(timeS)
    % Ensure gasTemperatureK is indexed correctly
    currentGasTempK = gasTemperatureK(i);
    if currentGasTempK > initialBarrelTempK && simulationResults.projectileVelocityMps(i) > 1e-3 % Only if gas is hot and projectile is moving
        currentPositionM = max(0, projectilePositionM(i)); % Ensure non-negative position for area calc
        exposedCylinderAreaM2 = pi * boreDiameterM * currentPositionM;
        contactAreaM2 = estimatedChamberAreaM2 + exposedCylinderAreaM2;
        heatLossRateW(i) = convectionCoefficient * contactAreaM2 * (currentGasTempK - initialBarrelTempK); % W
    end
end
heatLossRateMW = heatLossRateW / 1e6; % Rate in MW

% --- Create Figure ---
figureTitle = sprintf('Simulation Results: %s', parameters.loadedCaseName);
if isfield(parameters, 'propellantName') && ~isempty(parameters.propellantName)
    figureTitle = [figureTitle, sprintf(' (%s)', parameters.propellantName)];
end
figure('Name', figureTitle, 'NumberTitle', 'off', 'WindowState', 'maximized');

% --- Plot Layout (e.g., 4x3 or 5x2) ---
subplot(4,3,1); % Pressure vs Time
plot(timeMs, pressureMPa);
xlabel('Time [ms]'); ylabel('Pressure [MPa]');
title('Mean Gas Pressure'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,4); % Position vs Time
plot(timeMs, positionCm, 'DisplayName', 'Trajectory');
hold on;
plot(timeMs(end), positionCm(end), 'r*', 'MarkerSize', 8, 'DisplayName', 'End Sim.');
% --- Use CORRECTED parameter name ---
yline(parameters.barrelLength_m * 100, 'k--', 'DisplayName', 'Barrel Length');
% --- End CORRECTED parameter name ---
hold off;
xlabel('Time [ms]'); ylabel('Position [cm]');
title('Projectile Position'); legend('Location','northwest'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,7); % Linear Velocity vs Time
plot(timeMs, projectileVelocityMps);
xlabel('Time [ms]'); ylabel('Velocity [m/s]');
title('Linear Velocity'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,10); % Angular Velocity vs Time
plot(timeMs, omegaRpm);
xlabel('Time [ms]'); ylabel('Velocity [RPM]');
title('Angular Velocity'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,2); % Masses vs Time
plot(timeMs, propellantMassG, 'DisplayName', 'Solid');
hold on; plot(timeMs, gasMassG, '--', 'DisplayName', 'Gas'); hold off;
xlabel('Time [ms]'); ylabel('Mass [g]');
title('Propellant and Gas Masses'); legend('Location', 'best'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,5); % Temperature vs Time
plot(timeMs, gasTemperatureK);
xlabel('Time [ms]'); ylabel('Temperature [K]');
title('Mean Gas Temperature'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,8); % Bore Resistance Work vs Time
plot(timeMs, frictionWorkKJ);
xlabel('Time [ms]'); ylabel('Work [kJ]');
title('Cumulative Bore Resistance Work'); grid on; xlim([0, timeMs(end)]); % Title updated for clarity

subplot(4,3,11); % Cumulative Heat Loss vs Time
plot(timeMs, heatLossKJ);
xlabel('Time [ms]'); ylabel('Heat Loss [kJ]');
title('Cumulative Heat Loss to Barrel'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,3); % Torques vs Time (Rifling and Net)
plot(timeMs, riflingTorqueNm, 'DisplayName', 'T Rifling');
hold on;
plot(timeMs, netTorqueNm, '--', 'DisplayName', 'T Net');
hold off;
xlabel('Time [ms]'); ylabel('Torque [N*m]');
title('Main Torques'); legend('Location', 'best'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,6); % Torques vs Time (Resistive)
plot(timeMs, engravingTorqueNm, 'DisplayName', 'T Engraving');
hold on;
plot(timeMs, frictionTorqueNm, '--', 'DisplayName', 'T Rot. Friction');
plot(timeMs, engravingTorqueNm + frictionTorqueNm, ':', 'DisplayName', 'T Resist. Total');
hold off;
xlabel('Time [ms]'); ylabel('Torque [N*m]');
title('Resistive Torques'); legend('Location', 'best'); grid on; xlim([0, timeMs(end)]);

subplot(4,3,9); % Heat Loss Rate vs Time
plot(timeMs, heatLossRateMW);
xlabel('Time [ms]'); ylabel('Heat Loss Rate [MW]');
title('Instantaneous Heat Loss Rate'); grid on; xlim([0, timeMs(end)]);

% Empty plot to fill the last space (4,3,12)
subplot(4,3,12);
axis off; % Hide empty axes
% Add summary text here if desired
% text(0.1, 0.5, sprintf('V_{fin}: %.1f m/s\nRPM_{fin}: %.0f\nP_{max}: %.1f MPa', ...
%       simulationResults.projectileVelocityMps(end), omegaRpm(end), max(pressureMPa)), 'FontSize', 10);

% Overall Figure Title
sgtitle(figureTitle, 'Interpreter', 'none'); % Interpreter none to avoid issues with underscores

disp('Plots generated.');

end % End of function displaySimulationResults