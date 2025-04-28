% gui/displayEnergyWindow.m
% =========================================================
% Function to Display Energy Results, Pressure and BC in a Separate Window
% CORRECTED: Uses parameter names consistent with loadCaseParameters output.
% CORRECTED: Calculates energy balance vs. Released Powder Energy.
% MODIFIED: Uses normalized units and resize callback for proportional layout.
% MODIFIED: Increased space below Energy Balance title.
% MODIFIED: Adjusted label/value widths to prevent bold text clipping.
% MODIFIED: Split Energy Balance title into two lines for clarity.
% MODIFIED: Displays input G1 BC from parameters struct instead of calculated BC argument.
% MODIFIED: Added Muzzle Angular Velocity display.
% MODIFIED: Restored Estimated Barrel Temperature Increase display.
% =========================================================
function displayEnergyWindow(risultati_sim, parametri, delta_T_barrel, ~)
% INPUTS:
%   risultati_sim: Struct risultati simulazione (from runSimulation).
%   parametri:     Struct parametri simulazione (MUST contain 'bullet_bc_g1').
%   delta_T_barrel: Estimated barrel temp increase [K] (da calculateEnergyBalance).
%   ballistic_coefficient_calculated: Calculated BC value (e.g., G1) (da calculateEnergyBalance) - NOT used for BC display.

    disp('Creating/Updating responsive energy, pressure and BC results window...');

    % --- Input Validation (Basic) ---
    results_required = {'timeS', 'projectileVelocityMps', 'angularVelocityRadps', 'frictionWorkJ', 'heatLossJ', 'gasTemperatureK', 'remainingPropellantMassKg', 'gasPressurePa'};
    params_required = {'projMass_m', 'projMomentOfInertia_Ip', 'initialPropellantMass_m', 'specificGasConstant_R', 'specificHeatRatio_gamma', 'impetus_F', 'projWeight_gr', 'projDiameter_in', 'bullet_bc_g1'};
    valid_sim_data = isstruct(risultati_sim) && all(isfield(risultati_sim, results_required)) && ~isempty(risultati_sim.timeS);
    valid_params = isstruct(parametri) && all(isfield(parametri, params_required));

    if ~valid_sim_data || ~valid_params
        warning('Missing/invalid inputs for displayEnergyWindow. Check passed structs.');
         if ~valid_sim_data, disp('Error: risultati_sim struct is invalid or missing fields.'); disp('Required:'); disp(results_required'); if isstruct(risultati_sim), disp('Present:'); disp(fieldnames(risultati_sim)); end; end
         if ~valid_params, disp('Error: parametri struct is invalid or missing fields.'); disp('Required:'); disp(params_required'); if isstruct(parametri), disp('Present:'); disp(fieldnames(parametri)); end; end
        error('Cannot create energy window due to invalid input data.');
    end

    % --- Retrieve Input G1 BC for display ---
    % input_bc_g1 = parametri.bullet_bc_g1;
    % if ~isnumeric(input_bc_g1) || isnan(input_bc_g1)
    %      bc_display_string = 'N/A';
    % else
    %      bc_display_string = sprintf('%.3f', input_bc_g1);
    % end

    % --- Calculate Final Values and Energies ---
    v_muzzle=NaN; P_max_final_MPa=NaN; KE_linear=NaN; KE_rotational=NaN; Work_friction_total=NaN;
    Q_loss_total=NaN; U_gas_final=NaN; E_powder_released=NaN; m_burnt=NaN; perc_KE_lin=NaN;
    perc_KE_rot=NaN; perc_W_fric=NaN; perc_Q_loss=NaN; perc_U_gas=NaN; m_i=NaN; mass_burnt_perc_str = '';
    omega_muzzle = NaN; omega_final_rpm = NaN; % Added for omega
    try
        idx_end = length(risultati_sim.timeS);
        v_muzzle = risultati_sim.projectileVelocityMps(idx_end);
        omega_muzzle = risultati_sim.angularVelocityRadps(idx_end); % Get omega in rad/s
        omega_final_rpm = omega_muzzle * (60 / (2 * pi)); % Calculate RPM
        T_gas_final = risultati_sim.gasTemperatureK(idx_end);
        m_prop_rem_final = max(0, risultati_sim.remainingPropellantMassKg(idx_end));
        Work_friction_total = risultati_sim.frictionWorkJ(idx_end);
        Q_loss_total = risultati_sim.heatLossJ(idx_end);
        if isfield(risultati_sim, 'gasPressurePa') && ~isempty(risultati_sim.gasPressurePa) && all(~isnan(risultati_sim.gasPressurePa))
            P_max_final_MPa = max(risultati_sim.gasPressurePa) / 1e6;
        else
            warning('Max Pressure calculation skipped or invalid (NaN/empty).'); P_max_final_MPa = NaN;
        end
        m_p = parametri.projMass_m; I_p = parametri.projMomentOfInertia_Ip; m_i = parametri.initialPropellantMass_m;
        R_gas = parametri.specificGasConstant_R; gamma_i = parametri.specificHeatRatio_gamma; F_i = parametri.impetus_F;
        m_burnt = m_i - m_prop_rem_final;
        m_gas_final = max(1e-12, m_burnt);
        KE_linear = 0.5 * m_p * v_muzzle^2;
        KE_rotational = 0.5 * I_p * omega_muzzle^2;
        Cv=NaN; U_gas_final=NaN;
        if ~isnan(gamma_i) && (gamma_i - 1) > 1e-6 && ~isnan(R_gas) && R_gas > 0
            Cv = R_gas / (gamma_i - 1);
            if ~isnan(T_gas_final)
                 U_gas_final = m_gas_final * Cv * T_gas_final;
            end
        end
        E_powder_released = max(0, m_burnt * F_i);
        if ~isnan(m_burnt) && ~isnan(m_i) && m_i > 1e-9
             mass_burnt_perc_str = sprintf('[Mass Burned: %.2f%%]', (m_burnt/m_i)*100);
        end
        if ~isnan(E_powder_released) && E_powder_released > 1e-6
            if ~isnan(KE_linear), perc_KE_lin = (KE_linear / E_powder_released) * 100; end
            if ~isnan(KE_rotational), perc_KE_rot = (KE_rotational / E_powder_released) * 100; end
            if ~isnan(Work_friction_total), perc_W_fric = (Work_friction_total / E_powder_released) * 100; end
            if ~isnan(Q_loss_total), perc_Q_loss = (Q_loss_total / E_powder_released) * 100; end
            if ~isnan(U_gas_final), perc_U_gas = (U_gas_final / E_powder_released) * 100; end
        end
    catch ME_calc
        warning('Error calculating energy/pressure values inside displayEnergyWindow: %s. Some values may be NaN.', ME_calc.message);
    end

    % --- Check if Figure Exists, Create or Bring to Front ---
    figTag = 'EnergySummaryWindow';
    hFig = findobj('Type', 'figure', 'Tag', figTag);
    if isempty(hFig)
        disp('Creating new Energy Summary window.');
        % Adjusted height slightly to accommodate the new row
        figWidthPixels = 450; figHeightPixels = 395; % Increased height further for DeltaT
        screenSize = get(0, 'ScreenSize');
        figPosX = min( screenSize(3)-figWidthPixels, max(1, screenSize(3) * 0.6) );
        figPosY = (screenSize(4) - figHeightPixels) / 2;
        hFig = figure('Name', 'Energy, Pressure & BC Summary', ...
                      'Units', 'pixels', 'Position', [figPosX, figPosY, figWidthPixels, figHeightPixels], ...
                      'Units', 'normalized', ...
                      'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
                      'Resize', 'on', 'Color', [0.94 0.94 0.94], ...
                      'Tag', figTag, 'Visible', 'off', ...
                      'ResizeFcn', @energyWindowResizeFcn); % Attach resize function
        handles = struct(); handles.figure = hFig;
        % Create all text controls (ADDING Omega, RESTORING DeltaT controls)
        handles.titleMain = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Main Results', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11, 'Units', 'normalized','BackgroundColor',get(hFig,'Color'));
        handles.labelVel = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Muzzle Velocity:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueVel = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueVel');
        % --- ADDED Omega Controls ---
        handles.labelOmega = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Muzzle Angular Velocity:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); % New Label
        handles.valueOmega = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueOmega'); % New Value Field
        % --- END ADDED Omega Controls ---
        handles.labelPmax = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Maximum Pressure:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valuePmax = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valuePmax');
        % handles.labelBC = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Ballistic Coefficient (G1):', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueBC = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueBC');
        handles.titleEnergyBase = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Energy Balance', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11, 'Units', 'normalized','BackgroundColor',get(hFig,'Color'));
        handles.titleEnergySub = uicontrol('Parent', hFig, 'Style', 'text', 'String', '', 'HorizontalAlignment', 'center', 'FontWeight', 'normal', 'FontSize', 10, 'Units', 'normalized','BackgroundColor',get(hFig,'Color'));
        handles.labelKELin = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Linear Kinetic Energy:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueKELin = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueKELin');
        handles.labelKERot = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Rotational Kinetic Energy:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueKERot = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueKERot');
        handles.labelWorkBore = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Work Lost Bore Resistance:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueWorkBore = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueWorkBore');
        handles.labelQLoss = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Heat Lost to Barrel:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueQLoss = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueQLoss');
        handles.labelUGas = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Residual Internal Gas Energy:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color')); handles.valueUGas = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueUGas');
        % --- RESTORED DeltaT Controls ---
        handles.labelDeltaT = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Est. Barrel Temp Increase:', 'HorizontalAlignment', 'right', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'));
        handles.valueDeltaT = uicontrol('Parent', hFig, 'Style', 'text', 'String', 'N/A', 'HorizontalAlignment', 'left', 'Units', 'normalized','BackgroundColor',get(hFig,'Color'), 'Tag', 'valueDeltaT');
        % --- END RESTORED DeltaT Controls ---
        guidata(hFig, handles); % Save handles
    else
        disp('Found existing Energy Summary window. Updating...');
        figure(hFig); % Bring to front
        handles = guidata(hFig); % Retrieve handles
    end

    % --- Update Text Values ---
    formatValuePerc = @(value_J, perc) iff(~isnan(value_J) && ~isnan(perc), sprintf('%7.3f kJ (%5.1f %%)', value_J / 1000, perc), 'N/A');
    energyTitleBaseStr = 'Energy Balance';
    energyTitleSubStr = sprintf('(vs Released Energy ~ %.1f kJ)', E_powder_released / 1000);
    if isnan(E_powder_released), energyTitleSubStr = '(vs Released Energy N/A)'; end

    set(handles.titleEnergyBase, 'String', energyTitleBaseStr);
    set(handles.titleEnergySub, 'String', [energyTitleSubStr, ' ', mass_burnt_perc_str]);
    set(handles.valueVel, 'String', iff(isnan(v_muzzle), 'N/A', sprintf('%.1f m/s', v_muzzle)));
    % --- ADDED Omega Update ---
    set(handles.valueOmega, 'String', iff(isnan(omega_final_rpm), 'N/A', sprintf('%.0f RPM', omega_final_rpm)));
    % --- END ADDED Omega Update ---
    set(handles.valuePmax, 'String', iff(isnan(P_max_final_MPa), 'N/A', sprintf('%.2f MPa', P_max_final_MPa)));
    % set(handles.valueBC, 'String', bc_display_string); % Uses input BC G1
    set(handles.valueKELin, 'String', formatValuePerc(KE_linear, perc_KE_lin));
    set(handles.valueKERot, 'String', formatValuePerc(KE_rotational, perc_KE_rot));
    set(handles.valueWorkBore, 'String', formatValuePerc(Work_friction_total, perc_W_fric));
    set(handles.valueQLoss, 'String', formatValuePerc(Q_loss_total, perc_Q_loss));
    set(handles.valueUGas, 'String', formatValuePerc(U_gas_final, perc_U_gas));
    % --- RESTORED DeltaT Update ---
    set(handles.valueDeltaT, 'String', iff(isnan(delta_T_barrel), 'N/A', sprintf('%.1f K', delta_T_barrel)));
    % --- END RESTORED DeltaT Update ---


    % --- Apply Layout & Make Visible ---
    energyWindowResizeFcn(hFig, []); % Call resize function to position elements
    set(hFig, 'Visible', 'on'); % Make visible/bring to front
    disp('Energy summary window created/updated.');

end % End function displayEnergyWindow


% --- Resize Callback Function ---
function energyWindowResizeFcn(hObject, ~)
% Callback per ridimensionare e riposizionare i controlli.
% Adjusted to include Omega and DeltaT rows.

handles = guidata(hObject); if isempty(handles), return; end

% Define Layout Parameters (Normalized) - Adjust these if needed
marginX = 0.05; marginY = 0.025; % Reduced margin further
labelWidth = 0.52; valueWidth = 0.39; labelValueGap = 0.01;
itemHeight = 0.055; % Reduced height per item
titleHeight = 0.06; % Reduced title height
subtitleHeight = 0.045; % Reduced subtitle height
sectionGap = 0.025; % Reduced gap
currentY = 1.0 - marginY; % Start from top margin

% Helper functions for positioning
setpos = @(h, y, hgt, align) set(h, 'Position', [marginX, y - hgt, 1-2*marginX, hgt], 'HorizontalAlignment', align);
setpair = @(hL, hV, y, hgt) {set(hL, 'Position', [marginX, y - hgt, labelWidth, hgt]), ...
                           set(hV, 'Position', [marginX+labelWidth+labelValueGap, y - hgt, valueWidth, hgt]) };

try % Wrap layout updates in try-catch
    % --- Position Elements ---
    setpos(handles.titleMain, currentY, titleHeight, 'center');
    currentY = currentY - titleHeight - 0.01;

    setpair(handles.labelVel, handles.valueVel, currentY, itemHeight);
    currentY = currentY - itemHeight; % Move down for the next item

    % --- ADDED POSITIONING FOR OMEGA ---
    setpair(handles.labelOmega, handles.valueOmega, currentY, itemHeight); % Position Omega
    currentY = currentY - itemHeight; % Move down again
    % --- END ADDED SECTION ---

    setpair(handles.labelPmax, handles.valuePmax, currentY, itemHeight);
    currentY = currentY - itemHeight;

    % setpair(handles.labelBC, handles.valueBC, currentY, itemHeight);
    % currentY = currentY - itemHeight;

    currentY = currentY - sectionGap; % Gap after Main Results

    % --- Position Energy Balance section ---
    setpos(handles.titleEnergyBase, currentY, titleHeight, 'center');
    currentY = currentY - titleHeight;
    setpos(handles.titleEnergySub, currentY, subtitleHeight, 'center');
    currentY = currentY - subtitleHeight - 0.015; % Gap after Energy Title

    setpair(handles.labelKELin, handles.valueKELin, currentY, itemHeight); currentY = currentY - itemHeight;
    setpair(handles.labelKERot, handles.valueKERot, currentY, itemHeight); currentY = currentY - itemHeight;
    setpair(handles.labelWorkBore, handles.valueWorkBore, currentY, itemHeight); currentY = currentY - itemHeight;
    setpair(handles.labelQLoss, handles.valueQLoss, currentY, itemHeight); currentY = currentY - itemHeight;
    setpair(handles.labelUGas, handles.valueUGas, currentY, itemHeight);
    currentY = currentY - itemHeight; % Move down

    % --- RESTORED POSITIONING FOR DELTA T ---
    setpair(handles.labelDeltaT, handles.valueDeltaT, currentY, itemHeight); % Position Delta T at the end
    % --- END RESTORED POSITIONING ---

catch ME_resize
     warning('Error during energy window resize: %s', ME_resize.message);
     % Avoid erroring out completely during resize
end
end % End energyWindowResizeFcn


% --- Local helper function ---
function out = iff(condition, trueVal, falseVal)
    if condition, out = trueVal; else out = falseVal; end
end