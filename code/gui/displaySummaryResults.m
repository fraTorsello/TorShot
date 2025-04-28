function displaySummaryResults(hFigure)
% Displays a summary of the simulation results, including energy breakdown.
% CORRECTED: Uses parameter names consistent with loadCaseParameters output.
% Uses fprintf/disp for output.

    % --- Retrieve Data ---
    try
        handles = guidata(hFigure);
        if isempty(handles.risultati_sim) || ~isstruct(handles.risultati_sim) || ~isfield(handles.risultati_sim, 'timeS') || isempty(handles.risultati_sim.timeS) % Check timeS
            disp('[Helper] displaySummaryResults: Cannot show summary: no results.');
            return
        end
        res = handles.risultati_sim; % Uses runSimulation output names
        params = handles.parametri; % Uses loadCaseParameters output names
    catch ME_data
         disp('[Helper] displaySummaryResults: Error accessing handles data.');
         disp(ME_data.message);
         return
    end

    % --- Verify Required Fields ---
    % Results fields (from runSimulation output)
    required_res = {'timeS', 'projectileVelocityMps', 'angularVelocityRadps', 'gasPressurePa', ...
                    'remainingPropellantMassKg', 'frictionWorkJ', 'heatLossJ', 'gasTemperatureK', 'gasMassKg'};
    % Parameter fields (from loadCaseParameters output)
    required_params = {'initialPropellantMass_m', 'projMass_m', 'projMomentOfInertia_Ip', ...
                       'specificGasConstant_R', 'specificHeatRatio_gamma', 'impetus_F'}; % Using impetus_F for energy

    if ~all(isfield(res, required_res)) || ~all(isfield(params, required_params))
        disp('[Helper] displaySummaryResults: Cannot show full energy summary: missing data in results or parameters.');
        % Show basic summary if possible
        try
            idx_end = length(res.timeS);
            v_final = res.projectileVelocityMps(idx_end);
            omega_final_rpm = res.angularVelocityRadps(idx_end) * (60 / (2 * pi));
            p_max = max(res.gasPressurePa);
            t_exit = res.timeS(idx_end);
            prop_consumed_frac = 0;
            % --- Use CORRECT parameter name ---
            initialMassParam = params.initialPropellantMass_m;
            % --- End CORRECT parameter name ---
            if initialMassParam > 0
                 prop_consumed_frac = (initialMassParam - res.remainingPropellantMassKg(idx_end)) / initialMassParam;
            end

            fprintf('\n--- Simulation Results Summary (BASIC) ---\n');
            fprintf(' Final Velocity: %.2f m/s\n', v_final);
            fprintf(' Final RPM: %.0f RPM\n', omega_final_rpm);
            fprintf(' Max Pressure: %.2f MPa\n', p_max / 1e6);
            fprintf(' Exit Time: %.4f ms\n', t_exit * 1000);
            fprintf(' Propellant Consumed Fraction: %.2f %%\n', prop_consumed_frac * 100);
            if isfield(res, 'frictionWorkJ'), fprintf(' Total Bore Resistance Work: %.3f kJ\n', res.frictionWorkJ(idx_end) / 1000); end % Use frictionWorkJ
            if isfield(res, 'heatLossJ'), fprintf(' Total Heat Loss: %.3f kJ\n', res.heatLossJ(idx_end) / 1000); end % Use heatLossJ
            fprintf('-------------------------------------------\n');
        catch basicErr
             fprintf(2, '[Helper] displaySummaryResults: Error generating basic summary: %s\n', basicErr.message);
        end
        return; % Exit if full data missing
    end

    disp('[Helper] displaySummaryResults: Displaying full summary (in Command Window)...');

    % --- Base Calculations ---
    idx_end = length(res.timeS);
    t_exit = res.timeS(idx_end);
    p_max = max(res.gasPressurePa);

    % Final values from results
    v_final = res.projectileVelocityMps(idx_end);
    omega_final_rads = res.angularVelocityRadps(idx_end);
    m_prop_rem_final = res.remainingPropellantMassKg(idx_end);
    T_gas_final = res.gasTemperatureK(idx_end);
    W_fric_final = res.frictionWorkJ(idx_end); % Use frictionWorkJ (Bore Res Work)
    Q_loss_final = res.heatLossJ(idx_end);
    m_gas_final = res.gasMassKg(idx_end); % Use gasMassKg

    % --- Use CORRECT parameter names ---
    initialMass = params.initialPropellantMass_m;
    projMass = params.projMass_m;
    projInertia = params.projMomentOfInertia_Ip;
    gasConst = params.specificGasConstant_R;
    gamma = params.specificHeatRatio_gamma;
    impetus = params.impetus_F; % Use Impetus for energy
    % --- End CORRECT parameter names ---

    % Derived values
    omega_final_rpm = omega_final_rads * (60 / (2 * pi));
    prop_consumed_kg = 0; prop_consumed_frac = 0;
    if initialMass > 0
        prop_consumed_kg = initialMass - m_prop_rem_final;
        prop_consumed_frac = prop_consumed_kg / initialMass;
    end
    prop_consumed_kg = max(0, prop_consumed_kg); % Ensure non-negative

    % --- Energy Calculations ---
    KE_lin = 0.5 * projMass * v_final^2;         % J
    KE_rot = 0.5 * projInertia * omega_final_rads^2; % J

    Cv = 0; U_gas_final = 0;
    if (gamma - 1) > 1e-6 % Avoid division by zero
        Cv = gasConst / (gamma - 1); % J/(kg*K)
        U_gas_final = m_gas_final * Cv * T_gas_final; % Final gas internal energy [J]
    end

    % Powder Energy based on consumed mass and impetus
    E_powder_released = prop_consumed_kg * impetus; % [J]

    % --- Calculate Percentages ---
    perc_KE_lin = 0; perc_KE_rot = 0; perc_W_fric = 0;
    perc_Q_loss = 0; perc_U_gas = 0;

    if E_powder_released > 1e-6 % Avoid division by zero
        perc_KE_lin = (KE_lin / E_powder_released) * 100;
        perc_KE_rot = (KE_rot / E_powder_released) * 100;
        perc_W_fric = (W_fric_final / E_powder_released) * 100; % W_fric_final is Bore Res Work
        perc_Q_loss = (Q_loss_final / E_powder_released) * 100;
        perc_U_gas  = (U_gas_final / E_powder_released) * 100;
    end

    % --- Print Full Summary ---
    fprintf('\n--- Simulation Results Summary ---\n');
    fprintf(' Final Velocity:               %.2f m/s\n', v_final);
    fprintf(' Final RPM:                    %.0f RPM\n', omega_final_rpm);
    fprintf(' Max Pressure:                 %.2f MPa\n', p_max / 1e6);
    fprintf(' Exit Time:                    %.4f ms\n', t_exit * 1000);
    fprintf(' Propellant Consumed:        %.2f %%\n', prop_consumed_frac * 100);
    fprintf('\n--- Estimated Energy Balance (vs Released Energy = %.2f kJ) ---\n', E_powder_released / 1000);
    fprintf(' -> Linear Kinetic Energy:      %6.3f kJ (%5.1f %%)\n', KE_lin / 1000, perc_KE_lin);
    fprintf(' -> Rotational Kinetic Energy:  %6.3f kJ (%5.1f %%)\n', KE_rot / 1000, perc_KE_rot);
    fprintf(' -> Bore Resistance Work Loss:  %6.3f kJ (%5.1f %%)\n', W_fric_final / 1000, perc_W_fric); % Renamed W_fric
    fprintf(' -> Heat Loss to Barrel:        %6.3f kJ (%5.1f %%)\n', Q_loss_final / 1000, perc_Q_loss);
    fprintf(' -> Residual Gas Internal Energy:%6.3f kJ (%5.1f %%)\n', U_gas_final / 1000, perc_U_gas);
    % Note: Unaccounted energy includes residual gas KE, other friction/engraving not in W_fric_final, model approximations.
    fprintf('-----------------------------------------------------------\n');

end