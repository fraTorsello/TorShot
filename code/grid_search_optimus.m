% =========================================================
% Script: grid_search_optimus_v2.m
% VERSION: From scratch - Uses loadComponentData, external physics config,
%          pressure-dependent h_conv, direct result extraction.
% Descrizione: Esegue una ricerca su griglia per burnRateCoeff_a e burnRateExponent_beta
%              MINIMIZZANDO l'errore rispetto a una velocità target,
%              con P_max in range specificato.
% Metodo: Load Components -> Load Physics -> Assemble Params -> Grid Search -> Direct Analysis.
% =========================================================
clear; clc; close all;
fprintf('*** STARTING Grid Search Script (V2 - From Scratch) - %s ***\n', datestr(now));

% --- Constants and Path Setup ---
R_univ = 8.314462; % Universal gas constant [J/(mol*K)]
C = struct( ...
    'KG_PER_GRAIN', 1/15432.35835, ...
    'METERS_PER_INCH', 0.0254, ...
    'MM_PER_M', 1000, ...
    'G_PER_KG', 1000, ...
    'M3_PER_GRH2O', (1/15432.35835)/1000.0, ...
    'CM3_PER_M3', 1e6, ...
    'MM2_PER_M2', 1e6, ...
    'MPA_PER_PA', 1e-6, ...
    'PA_PER_MPA', 1e6, ...
    'GRAINS_PER_KG', 15432.35835, ...
    'INCHES_PER_METER', 1/0.0254, ...
    'GRH2O_PER_M3', 1000.0*15432.35835, ...
    'IN2_PER_M2', (1/0.0254)^2 ...
);

disp("Adding paths...");
try
    script_path = fileparts(mfilename('fullpath'));
    if isempty(script_path), script_path = pwd; end
    addpath(fullfile(script_path));
    addpath(genpath(fullfile(script_path, 'config')));
    addpath(genpath(fullfile(script_path, 'core')));
    addpath(genpath(fullfile(script_path, 'utils')));
    addpath(genpath(fullfile(script_path, 'postprocessing')));
    addpath(genpath(fullfile(script_path, 'gui')));
    rehash path;
    disp("Paths OK.");
    % Check for required functions
    if exist('loadComponentData', 'file') ~= 2, error('loadComponentData.m not found.'); end
    if exist('runSimulation', 'file') ~= 2, error('runSimulation.m not found.'); end
    % Define getfield_safe locally if potentially missing from path
    if exist('getfield_safe', 'file') ~= 2
         warning('getfield_safe.m not found on path, defining locally.');
         getfield_safe = @(s, field, default) local_getfield_safe(s, field, default);
    end
catch ME_path, warning('Path error: %s', ME_path.message); return; end

% --- Define Base Case Components and Inputs ---
% <<< MODIFY THESE BASED ON THE DESIRED BASE CASE >>>
baseCaseDescription = '.243 Win / Varget / 95gr Nosler Partition';
cartridgeSelection = 'cart_243_Winchester';
bulletSelection    = 'cal_243_95gr_Nosler_Partition_SP';
powderSelection    = 'varget_params';

% --- User Inputs for Base Case (Mimic GUI Input) ---
inputs = struct();
% !!! VERIFY THESE VALUES MATCH YOUR INTENDED BASE CASE !!!
inputs.PropellantCharge_gr = 39;  
inputs.SeatingDepth_in     = 0.365; % Reference value, ensure matches bullet file if used
inputs.BarrelLength_in     = 24.0;
inputs.TwistRate_inPerTurn = 12.0;
inputs.InitialTemperature_K= 293.15;% Standard temp (20 C)
inputs.AmbientPressure_Pa  = 101325;% Standard pressure
inputs.ShotStartPressure_MPa = 25.0; % Standard estimate

fprintf('--- Loading Base Case Components for: %s ---\n', baseCaseDescription);
fprintf('  Cartridge: %s\n', cartridgeSelection);
fprintf('  Bullet:    %s\n', bulletSelection);
fprintf('  Powder:    %s\n', powderSelection);
fprintf('  Using Base Inputs: Charge=%.1f gr, Barrel=%.1f in, Twist=%.1f\n', ...
        inputs.PropellantCharge_gr, inputs.BarrelLength_in, inputs.TwistRate_inPerTurn);

try
    % --- Load Component Data ---
    componentData = loadComponentData(cartridgeSelection, bulletSelection, powderSelection);
    cd = componentData;

    % --- Load Physics Configuration ---
    fprintf('Loading physics configuration...\n');
    physicsConfigFile = fullfile('config', 'defaultPhysicsConfig.m'); % Adjust path if needed
    if ~exist(physicsConfigFile, 'file'), error('Physics configuration file not found: %s', physicsConfigFile); end
    try
        clear physicsConfig; run(physicsConfigFile);
        if ~exist('physicsConfig', 'var') || ~isstruct(physicsConfig), error('Script "%s" did not define "physicsConfig" struct.', physicsConfigFile); end
        fprintf('Loaded physics configuration from %s\n', physicsConfigFile);
    catch ME_physicsLoad, error('Error loading/running physics config file "%s": %s', physicsConfigFile, ME_physicsLoad.message); end

    % --- Assemble Base Parameters ---
    fprintf('Assembling base simulation parameters...\n');
    % Convert GUI-like inputs to SI
    propCharge_kg = inputs.PropellantCharge_gr * C.KG_PER_GRAIN;
    barrelLength_m = inputs.BarrelLength_in * C.METERS_PER_INCH;
    twistRate_inPerTurn = inputs.TwistRate_inPerTurn;
    initialTemp_K = inputs.InitialTemperature_K;
    ambientPress_Pa = inputs.AmbientPressure_Pa;
    shotStartPressure_Pa = inputs.ShotStartPressure_MPa * C.PA_PER_MPA;

    % Get SI values directly from loaded component data
    bulletMass_kg = getfield_safe(cd.bullet, 'mass_kg', NaN);
    bulletDiameter_m = getfield_safe(cd.bullet, 'diameter_m', NaN);
    bulletLength_m = getfield_safe(cd.bullet, 'length_m', NaN);
    seatingDepth_m = getfield_safe(cd.bullet, 'seatingDepth_m', NaN); % Use loaded seating depth
    caseLength_m = getfield_safe(cd.cartridge, 'caseLength_m', NaN);
    maxCapacity_m3 = getfield_safe(cd.cartridge, 'maxCaseCapacity_m3', NaN);
    boreDiameter_m = getfield_safe(cd.cartridge, 'boreDiameter_m', NaN);
    powderDensity_kgm3 = getfield_safe(cd.powder, 'propellantDensity_rho_s', NaN);

    % Validate essential loaded SI values
    if isnan(bulletMass_kg) || isnan(bulletDiameter_m) || isnan(bulletLength_m) || isnan(seatingDepth_m) || ...
       isnan(caseLength_m) || isnan(maxCapacity_m3) || isnan(boreDiameter_m) || isnan(powderDensity_kgm3) || powderDensity_kgm3 <= 0
        error('One or more required values missing or invalid from loaded component data.');
    end
    fprintf('  Loaded Seating Depth = %.3f in (%.4f m)\n', seatingDepth_m * C.INCHES_PER_METER, seatingDepth_m);

    % Calculate Volumes
    powderVolume_m3 = 0; if propCharge_kg > 0, powderVolume_m3 = propCharge_kg / powderDensity_kgm3; end
    bulletArea_m2 = pi * (bulletDiameter_m / 2)^2;
    seatedBulletVolume_m3 = bulletArea_m2 * seatingDepth_m;
    qlUseableCapacity_m3 = maxCapacity_m3 - seatedBulletVolume_m3;
    initialFreeVolumeV0_m3 = qlUseableCapacity_m3 - powderVolume_m3;
    if qlUseableCapacity_m3 <= 0, warning('Calculated Useable Case Capacity <= 0.'); end
    if initialFreeVolumeV0_m3 <= 1e-12 && propCharge_kg > 0, error('Calculated Initial Free Volume V0 <= 0.'); end
    if propCharge_kg == 0, initialFreeVolumeV0_m3 = qlUseableCapacity_m3; end

    % Assemble the final baseParameters struct
    baseParameters = struct();
    % Identifiers
    baseParameters.loadedCartridgeName = getfield_safe(cd.cartridge, 'cartridgeName', cartridgeSelection);
    baseParameters.loadedBulletName = getfield_safe(cd.bullet, 'bulletName', bulletSelection);
    baseParameters.loadedPowderName = getfield_safe(cd.powder, 'powderName', powderSelection);
    % Original Units for reference
    baseParameters.projWeight_gr = bulletMass_kg * C.GRAINS_PER_KG;
    baseParameters.projDiameter_in = bulletDiameter_m * C.INCHES_PER_METER;
    % Core SI Parameters
    baseParameters.projMass_m = bulletMass_kg;
    baseParameters.projArea_Ab = bulletArea_m2;
    baseParameters.projMomentOfInertia_Ip = 0.5 * bulletMass_kg * (bulletDiameter_m / 2)^2; % Estimate
    baseParameters.initialPropellantMass_m = propCharge_kg;
    baseParameters.initialFreeVolume_V0 = initialFreeVolumeV0_m3;
    baseParameters.barrelLength_m = barrelLength_m;
    baseParameters.boreDiameter_Db = boreDiameter_m;
    baseParameters.twistRate_inPerTurn = twistRate_inPerTurn;
    if twistRate_inPerTurn ~= 0
        meters_per_turn = abs(twistRate_inPerTurn) * C.METERS_PER_INCH;
        baseParameters.twistRate_rad_m = sign(twistRate_inPerTurn) * (2 * pi) / meters_per_turn;
    else
        baseParameters.twistRate_rad_m = 0;
    end

    % Copy powder fields
    powderDataLoaded = cd.powder; powderFields = fieldnames(powderDataLoaded);
    fprintf('Copying parameters from loaded powder data...\n');
    for k=1:length(powderFields), baseParameters.(powderFields{k}) = powderDataLoaded.(powderFields{k}); end

    % Calculate specificGasConstant_R
    if isfield(baseParameters, 'molarMass_M_gas') && isnumeric(baseParameters.molarMass_M_gas) && baseParameters.molarMass_M_gas > 0
        baseParameters.specificGasConstant_R = R_univ / baseParameters.molarMass_M_gas;
        fprintf('Calculated specificGasConstant_R = %.2f J/(kg*K)\n', baseParameters.specificGasConstant_R);
    else, error('Molar mass (molarMass_M_gas) missing or invalid.'); end

    % Calculate Initial Surface Area S0
    sigma_powder = getfield_safe(baseParameters,'specificSurfaceArea_sigma', NaN); % Use value already copied
    if isnan(sigma_powder) || sigma_powder <= 0, error('Powder specific surface area (sigma) invalid.'); end
    baseParameters.initialSurfaceArea_S0 = baseParameters.initialPropellantMass_m * sigma_powder; % Use sigma_powder here
    if baseParameters.initialPropellantMass_m == 0, baseParameters.initialSurfaceArea_S0 = 0;
    elseif baseParameters.initialSurfaceArea_S0 <= 0, baseParameters.initialSurfaceArea_S0 = 1e-9; warning('Calculated initialSurfaceArea_S0 <= 0, set to 1e-9.'); end

    % Copy G1 BC
    baseParameters.bullet_bc_g1 = getfield_safe(cd.bullet, 'bc_g1', NaN);
    baseParameters.formFactor_G1 = getfield_safe(cd.bullet, 'formFactor_g1', NaN);
    fprintf('Using bullet_bc_g1 = %.3f, formFactor_G1 = %.4f from loaded bullet data.\n', baseParameters.bullet_bc_g1, baseParameters.formFactor_G1);

    % Assign Physics Params from Loaded Config
    fprintf('Assigning physics parameters from loaded configuration...\n');
    baseParameters.engravingTorque_Nm = getfield_safe(physicsConfig, 'engravingTorque_Nm', NaN);
    baseParameters.engravingEndPosition_m = getfield_safe(physicsConfig, 'engravingEndPosition_m', NaN);
    baseParameters.rotationalFrictionCoeff_muRot = getfield_safe(physicsConfig, 'rotationalFrictionCoeff_muRot', NaN);
    baseParameters.h_conv_base_coefficient = getfield_safe(physicsConfig, 'h_conv_base_coefficient', NaN);
    baseParameters.h_conv_pressure_reference_Pa = getfield_safe(physicsConfig, 'h_conv_pressure_reference_Pa', NaN);
    baseParameters.h_conv_pressure_exponent = getfield_safe(physicsConfig, 'h_conv_pressure_exponent', NaN);
    baseParameters.maxSafetyTime_s = getfield_safe(physicsConfig, 'maxSafetyTime_s', 0.004);
    baseParameters.boreResistance_travel_m = getfield_safe(physicsConfig, 'boreResistance_travel_m', []);
    baseParameters.boreResistance_pressure_Pa = getfield_safe(physicsConfig, 'boreResistance_pressure_Pa', []);

    % Validate loaded physics parameters
    if isnan(baseParameters.engravingTorque_Nm) || isnan(baseParameters.engravingEndPosition_m) || ...
       isnan(baseParameters.rotationalFrictionCoeff_muRot) || ...
       isnan(baseParameters.h_conv_base_coefficient) || isnan(baseParameters.h_conv_pressure_reference_Pa) || isnan(baseParameters.h_conv_pressure_exponent) || ...
       baseParameters.h_conv_pressure_reference_Pa <= 0 || ...
       isempty(baseParameters.boreResistance_travel_m) || isempty(baseParameters.boreResistance_pressure_Pa) || ...
       length(baseParameters.boreResistance_travel_m) ~= length(baseParameters.boreResistance_pressure_Pa)
         error('Failed to load valid physics parameters from config file. Check defaultPhysicsConfig.m.');
    end

    % Add initial conditions and environmental params
    baseParameters.initialTemperature_K = initialTemp_K;
    baseParameters.ambientPressure = ambientPress_Pa;
    baseParameters.shotStartPressure_Pa = shotStartPressure_Pa;

    % --- Final Validation (Optional but recommended) ---
     required_for_odes = { ...
            'initialPropellantMass_m', 'initialTemperature_K', 'projMass_m', 'projArea_Ab', 'initialFreeVolume_V0', ...
            'propellantDensity_rho_s', 'covolume_b', 'impetus_F', ...
            'specificHeatRatio_gamma', 'specificGasConstant_R', ...
            'burnRateCoeff_a', 'burnRateExponent_beta', 'initialSurfaceArea_S0', 'formFunctionParam_theta', ...
            'projMomentOfInertia_Ip', 'twistRate_rad_m', ...
            'shotStartPressure_Pa', ...
            'boreResistance_travel_m', 'boreResistance_pressure_Pa', ...
            'engravingTorque_Nm', 'engravingEndPosition_m', ...
            'rotationalFrictionCoeff_muRot', ...
            'h_conv_base_coefficient', 'h_conv_pressure_reference_Pa', 'h_conv_pressure_exponent', ... % New h_conv params
            'boreDiameter_Db', ...
            'barrelLength_m', 'maxSafetyTime_s', 'ambientPressure' ...
            };
     missing_odes_fields = required_for_odes(~isfield(baseParameters, required_for_odes));
     if ~isempty(missing_odes_fields), error('Assembled baseParameters missing fields for simulationOdes: %s', strjoin(missing_odes_fields,', ')); end
     fprintf('Base parameters assembled and validated successfully.\n');

catch ME_load_assemble
    fprintf(2, 'Error loading components or assembling base parameters: %s\n', ME_load_assemble.message);
    if ~isempty(ME_load_assemble.stack), fprintf(2, 'File: %s, Line: %d\n', ME_load_assemble.stack(1).file, ME_load_assemble.stack(1).line); end
    rethrow(ME_load_assemble);
end

% --- Define Parameter Ranges and Grid Resolution ---
param_names = {'burnRateCoeff_a', 'burnRateExponent_beta'};
param1_field = param_names{1};
param2_field = param_names{2};
if ~isfield(baseParameters, param1_field) || ~isfield(baseParameters, param2_field), error('Search parameters not found in baseParameters.'); end
initial_param1 = baseParameters.(param1_field); initial_param2 = baseParameters.(param2_field);
fprintf('Grid Center Values: %s=%.4e, %s=%.4f\n', param1_field, initial_param1, param2_field, initial_param2);

percent_range = 0.1; % Range +/- 0.8%
delta_param1 = percent_range * initial_param1; delta_param2 = percent_range * initial_param2;
lb = [initial_param1 - delta_param1; initial_param2 - delta_param2]; ub = [initial_param1 + delta_param1; initial_param2 + delta_param2];
if lb(1) <= 0, lb(1) = eps; end % Ensure a > 0
beta_max = 1.0; if ub(2) > beta_max, ub(2) = beta_max; end; beta_min = 0.5; if lb(2) < beta_min, lb(2) = beta_min; end % Limit beta range
fprintf('Parameter search range: %s=[%.3e, %.3e], %s=[%.4f, %.4f]\n', param1_field, lb(1), ub(1), param2_field, lb(2), ub(2));

N1 = 45; N2 = 45; % Grid resolution (adjust as needed)
total_simulations_grid = N1 * N2;
fprintf('Grid resolution: %dx%d = %d simulations.\n', N1, N2, total_simulations_grid);
param1_values = linspace(lb(1), ub(1), N1); param2_values = linspace(lb(2), ub(2), N2);

% --- Define Pressure Constraints and Target Velocity ---
P_min_limit_MPa = 380; P_max_limit_MPa = 395.0; % Pressure Target
target_velocity_mps = 915; % Speed Target
fprintf('Pressure Constraints: [%.1f, %.1f] MPa\n', P_min_limit_MPa, P_max_limit_MPa);
fprintf('Target Velocity: %.1f m/s\n', target_velocity_mps);

% --- Initialize Results Table ---
fprintf('\n--- Initializing Results Table ---\n');
results_table = table('Size',[total_simulations_grid 6], ... % Use total_simulations_grid size
                      'VariableTypes',{'double','double','double','double','double','double'}, ...
                      'VariableNames',{param1_field, param2_field,'P_max_MPa','Combustion_Perc','V_final_mps', 'Velocity_Error'});
results_idx = 0; % Index for storing results
best_velocity_error = Inf; best_param1 = NaN; best_param2 = NaN; best_P_max = NaN; best_V_final = NaN; best_combustion = NaN;
feasible_points_count = 0;

% --- Execute Simulations on the Grid ---
fprintf('\n--- Starting Grid Simulations (%s) ---\n', datestr(now));
start_time = tic;
for i = 1:N1
    current_param1 = param1_values(i);
    for j = 1:N2
        current_param2 = param2_values(j);
        current_sim_index = (i-1)*N2 + j;

        % Update parameters for current grid point
        currentParameters = baseParameters;
        currentParameters.(param1_field) = current_param1;
        currentParameters.(param2_field) = current_param2;

        % Progress Logging
        if mod(current_sim_index, 20) == 1 || current_sim_index == total_simulations_grid
           elapsed_time = toc(start_time);
           fprintf('Sim G %d/%d (%s=%.3e,%s=%.4f) T:%.1fs\n', current_sim_index, total_simulations_grid, param1_field, current_param1, param2_field, current_param2, elapsed_time);
        end

        % --- Run Simulation and Extract Results Directly ---
        P_max_MPa = NaN; V_final_mps = NaN; velocity_error_current = NaN; perc_comb_current = NaN; % Reset
        results = []; % Clear previous results

        try
            % Run the simulation
            results = runSimulation(currentParameters);

            % Check if runSimulation returned valid results
            if isempty(results) || ~isstruct(results) || ~isfield(results, 'timeS') || isempty(results.timeS) || length(results.timeS) < 2 || ~isfield(results, 'gasPressurePa') || ~isfield(results,'projectileVelocityMps') || ~isfield(results,'remainingPropellantMassKg')
                 fprintf(2, 'WARN: runSimulation returned invalid/empty results for params: a=%.3e, beta=%.4f\n', current_param1, current_param2);
            else
                % Extract results directly (avoid analyzeAndSaveResults)
                pres_raw = results.gasPressurePa;
                vel_raw = results.projectileVelocityMps;
                prop_rem_raw = results.remainingPropellantMassKg;

                % Check for NaN/Inf in results needed
                if any(isnan(pres_raw)) || any(isinf(pres_raw)) || any(isnan(vel_raw)) || any(isinf(vel_raw)) || any(isnan(prop_rem_raw)) || any(isinf(prop_rem_raw))
                    fprintf(2, 'WARN: NaN/Inf found in simulation output arrays for params: a=%.3e, beta=%.4f\n', current_param1, current_param2);
                else
                    % Calculate metrics if results are numeric
                    P_max_MPa = max(pres_raw) / 1e6;
                    V_final_mps = vel_raw(end);
                    prop_consumed_kg = max(0, currentParameters.initialPropellantMass_m - prop_rem_raw(end));
                    if currentParameters.initialPropellantMass_m > 1e-12
                        perc_comb_current = (prop_consumed_kg / currentParameters.initialPropellantMass_m) * 100;
                    else
                        perc_comb_current = 0; % Avoid division by zero if initial mass is zero
                    end
                    velocity_error_current = abs(V_final_mps - target_velocity_mps);
                end % End check for NaN/Inf in results
            end % End check for valid simulation results struct
        catch ME_sim_loop
             fprintf(2, 'ERROR during simulation run for params a=%.3e, beta=%.4f: %s\n', current_param1, current_param2, ME_sim_loop.message);
             % Optionally display stack trace for debugging:
             % disp(ME_sim_loop.getReport('basic'));
        end % End try-catch for simulation run

        % --- Save results to table (even if NaN) ---
        results_idx = results_idx + 1; % Increment index
        if results_idx > height(results_table)
             warning('Results index exceeds table size. Check table initialization.'); % Safety check
        else
            results_table(results_idx,:) = {current_param1, current_param2, P_max_MPa, perc_comb_current, V_final_mps, velocity_error_current};
        end

        % --- Check Constraints and Update Best Result (only if results are valid numbers) ---
        if ~isnan(P_max_MPa) && ~isnan(V_final_mps) && ~isnan(velocity_error_current) && ~isnan(perc_comb_current)
            if P_max_MPa >= P_min_limit_MPa && P_max_MPa <= P_max_limit_MPa % Check pressure constraint
                feasible_points_count = feasible_points_count + 1;
                if velocity_error_current < best_velocity_error % Check if better than current best
                    fprintf('----> New Best! (Sim G %d) Verr:%.2f|V:%.2f|P:%.2f|Comb:%.1f%%|%s:%.3e|%s:%.4f\n', ...
                           current_sim_index, velocity_error_current, V_final_mps, P_max_MPa, perc_comb_current, param1_field, current_param1, param2_field, current_param2);
                    % Update best values
                    best_velocity_error = velocity_error_current;
                    best_param1 = current_param1; best_param2 = current_param2;
                    best_P_max = P_max_MPa; best_V_final = V_final_mps;
                    best_combustion = perc_comb_current;
                end
            end % End pressure constraint check
        end % End if valid results check

    end % End loop param2 (j)
end % End loop param1 (i)

results_table = results_table(1:results_idx,:); % Trim potentially unused rows if loop exited early

% --- Report Final Results ---
fprintf('\n--- Grid Search Finished (%s) ---\n', datestr(now));
total_time = toc(start_time);
fprintf('Total Execution Time (grid only): %.1f minutes\n', total_time/60);
fprintf('Simulations Run on Grid: %d\n', results_idx);
fprintf('Feasible Points Found (Pmax in [%.1f, %.1f] MPa): %d\n', P_min_limit_MPa, P_max_limit_MPa, feasible_points_count);

if feasible_points_count > 0 && ~isnan(best_param1)
    fprintf('\n=== BEST OVERALL RESULT FOUND (Target Velocity=%.1f m/s) ===\n', target_velocity_mps);
    fprintf('Optimal Parameters (%s / %s):\n', param1_field, param2_field);
    fprintf('  %s = %.6e\n', param1_field, best_param1);
    fprintf('  %s = %.6f\n', param2_field, best_param2);
    fprintf('Corresponding Results:\n');
    fprintf('  Final Velocity    = %.2f m/s (Error: %.2f m/s)\n', best_V_final, best_velocity_error);
    fprintf('  Maximum Pressure  = %.4f MPa (Constraint: [%.1f, %.1f] MPa)\n', best_P_max, P_min_limit_MPa, P_max_limit_MPa);
    fprintf('  Combustion        = %.4f %%\n', best_combustion);

    % --- Save final results and table ---
    try
        save_filename = sprintf('grid_search_v2_results_%s_%s_%s_%dx%d_%s.mat', ...
                                strrep(cartridgeSelection,'_QL',''), strrep(bulletSelection,'_','-'), strrep(powderSelection,'_params',''), ...
                                N1, N2, datestr(now,'yyyymmdd_HHMMSS'));
        save(save_filename, 'best_param1', 'best_param2', 'best_V_final', 'best_velocity_error', 'target_velocity_mps', ...
               'best_P_max', 'best_combustion', 'param1_values', 'param2_values', 'lb', 'ub', ...
               'P_min_limit_MPa', 'P_max_limit_MPa', 'results_table', 'total_time', ...
               'baseCaseDescription', 'cartridgeSelection', 'bulletSelection', 'powderSelection', 'inputs', ...
               'initial_param1', 'initial_param2'); % Removed center point initial results as they are in table row 1 now
        fprintf('Best results and full table saved in: %s\n', save_filename);
    catch ME_save, warning('Error saving final results: %s', ME_save.message); end

    % --- Visualization ---
    try
        fprintf('Launching visualization...\n');
        figure; hFig = gcf; hFig.Name = sprintf('Grid Search V2 Results (%dx%d) - Target %.1f m/s', N1, N2, target_velocity_mps);
        results_to_plot = results_table(~isnan(results_table.P_max_MPa) & ~isnan(results_table.V_final_mps),:); % Plot only valid points
        if isempty(results_to_plot), warning('No valid simulation results to plot.'); close(hFig); else
            fprintf('Plotting: %d valid simulation points.\n', height(results_to_plot));
            feasible_mask = results_to_plot.P_max_MPa >= P_min_limit_MPa & results_to_plot.P_max_MPa <= P_max_limit_MPa;
            infeasible_mask = ~feasible_mask; num_feasible = sum(feasible_mask); num_infeasible = sum(infeasible_mask);
            fprintf('Plotting: %d feasible, %d infeasible.\n', num_feasible, num_infeasible);
            param1_data = results_to_plot.(param1_field); param2_data = results_to_plot.(param2_field);
            best_p1_plot = best_param1; best_p2_plot = best_param2;

            % Plotting code (same as before, ensure best point markers only plot if best point exists)
             % Subplot 1: Velocity Error
            subplot(1, 3, 1); hold on; h_plots_verr = gobjects(0); plot_labels_verr = {};
            if num_infeasible > 0, h_infeas_verr = scatter(param1_data(infeasible_mask), param2_data(infeasible_mask), 20, results_to_plot.Velocity_Error(infeasible_mask), 'filled', 'MarkerFaceAlpha', 0.3); h_plots_verr(end+1) = h_infeas_verr; plot_labels_verr{end+1} = sprintf('Infeasible (%d)', num_infeasible); end
            if num_feasible > 0, h_feas_verr = scatter(param1_data(feasible_mask), param2_data(feasible_mask), 30, results_to_plot.Velocity_Error(feasible_mask), 'filled'); h_plots_verr(end+1) = h_feas_verr; plot_labels_verr{end+1} = sprintf('Feasible (%d)', num_feasible); end
            if ~isnan(best_p1_plot) && ~isnan(best_p2_plot), h_best_verr = scatter(best_p1_plot, best_p2_plot, 100, 'r', 'x', 'LineWidth', 2); h_plots_verr(end+1) = h_best_verr; plot_labels_verr{end+1} = 'Best Overall'; end
            hold off; colorbar; xlabel(param1_field,'Interpreter','none'); ylabel(param2_field,'Interpreter','none'); title(sprintf('Velocity Error (m/s) - Min=%.2f', best_velocity_error)); if ~isempty(h_plots_verr), legend(h_plots_verr, plot_labels_verr, 'Location','best','AutoUpdate','off'); end; grid on;
            % Subplot 2: Final Velocity
            subplot(1, 3, 2); hold on; h_plots_vfin = gobjects(0); plot_labels_vfin = {};
            if num_infeasible > 0, h_infeas_vfin = scatter(param1_data(infeasible_mask), param2_data(infeasible_mask), 20, results_to_plot.V_final_mps(infeasible_mask), 'filled', 'MarkerFaceAlpha', 0.3); h_plots_vfin(end+1) = h_infeas_vfin; plot_labels_vfin{end+1} = sprintf('Infeasible (%d)', num_infeasible); end
            if num_feasible > 0, h_feas_vfin = scatter(param1_data(feasible_mask), param2_data(feasible_mask), 30, results_to_plot.V_final_mps(feasible_mask), 'filled'); h_plots_vfin(end+1) = h_feas_vfin; plot_labels_vfin{end+1} = sprintf('Feasible (%d)', num_feasible); end
             if ~isnan(best_p1_plot) && ~isnan(best_p2_plot), h_best_vfin = scatter(best_p1_plot, best_p2_plot, 100, 'r', 'x', 'LineWidth', 2); h_plots_vfin(end+1) = h_best_vfin; plot_labels_vfin{end+1} = 'Best Overall'; end
            hold off; colorbar; xlabel(param1_field,'Interpreter','none'); ylabel(param2_field,'Interpreter','none'); title(sprintf('Final Velocity (m/s) - Best=%.1f', best_V_final)); if ~isempty(h_plots_vfin), legend(h_plots_vfin, plot_labels_vfin, 'Location','best','AutoUpdate','off'); end; grid on;
            % Subplot 3: Max Pressure
            subplot(1, 3, 3); hold on; h_plots_pmax = gobjects(0); plot_labels_pmax = {};
            if num_infeasible > 0, h_infeas_pmax = scatter(param1_data(infeasible_mask), param2_data(infeasible_mask), 20, results_to_plot.P_max_MPa(infeasible_mask), 'filled', 'MarkerFaceAlpha', 0.3); h_plots_pmax(end+1) = h_infeas_pmax; plot_labels_pmax{end+1} = sprintf('Infeasible (%d)', num_infeasible); end
            if num_feasible > 0, h_feas_pmax = scatter(param1_data(feasible_mask), param2_data(feasible_mask), 30, results_to_plot.P_max_MPa(feasible_mask), 'filled'); h_plots_pmax(end+1) = h_feas_pmax; plot_labels_pmax{end+1} = sprintf('Feasible (%d)', num_feasible); end
            if ~isnan(best_p1_plot) && ~isnan(best_p2_plot), h_best_pmax = scatter(best_p1_plot, best_p2_plot, 100, 'r', 'x', 'LineWidth', 2); h_plots_pmax(end+1) = h_best_pmax; plot_labels_pmax{end+1} = 'Best Overall'; end
            hold off; colorbar; xlabel(param1_field,'Interpreter','none'); ylabel(param2_field,'Interpreter','none'); title(sprintf('Max Pressure (MPa) - Best=%.1f MPa', best_P_max)); if ~isempty(h_plots_pmax), legend(h_plots_pmax, plot_labels_pmax, 'Location','best','AutoUpdate','off'); end; grid on;

            sgtitle(sprintf('Grid Search V2 (%dx%d) %s (%s/%s/%s) - Target %.1fm/s - Pmax[%.0f,%.0f]MPa', N1, N2, baseCaseDescription, cartridgeSelection, bulletSelection, powderSelection, target_velocity_mps, P_min_limit_MPa, P_max_limit_MPa), 'Interpreter', 'none');
        end
    catch ME_plot, warning('Error during plot visualization: %s', ME_plot.message); if exist('hFig', 'var') && ishandle(hFig), try close(hFig); catch; end; end; end
else
     fprintf('\n!!! WARNING: No FEASIBLE points found (Pmax in range [%.1f, %.1f] MPa).\n', P_min_limit_MPa, P_max_limit_MPa);
     try % Save table anyway
        save_filename_nores = sprintf('grid_search_v2_NO_RESULTS_%s_%s_%s_%dx%d_%s.mat', ...
                                      strrep(cartridgeSelection,'_QL',''), strrep(bulletSelection,'_','-'), strrep(powderSelection,'_params',''), N1, N2, datestr(now,'yyyymmdd_HHMMSS'));
        save(save_filename_nores, 'param1_values', 'param2_values', 'lb', 'ub', 'P_min_limit_MPa', 'P_max_limit_MPa', 'target_velocity_mps', ...
               'results_table', 'total_time', 'baseCaseDescription', 'cartridgeSelection', 'bulletSelection', 'powderSelection', 'inputs'); % Removed center point initial results
        fprintf('Full simulation table (no feasible results) saved in: %s\n', save_filename_nores);
     catch ME_save_nores, warning('Error saving no-results table: %s', ME_save_nores.message); end
end

disp("===================================================");
disp("Grid Search Script (V2 - From Scratch) finished.");

% --- Local helper function definitions ---
function v = local_getfield_safe(s, field, default)
% Returns field value or default if field doesn't exist or is empty.
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end

% This might not be strictly necessary if getfield_safe is defined elsewhere or on path
function v = iffHelper(condition, trueVal, falseVal)
% Simple inline-if function substitute.
    if condition, v = trueVal; else, v = falseVal; end
end