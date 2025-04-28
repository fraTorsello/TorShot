% gui/plotResults.m
function plotResults(hFigure)
% Plots the simulation results stored in handles.risultati_sim onto the GUI axes.
% MODIFIED: Added subplots for barrel stresses and safety factor.

handles = guidata(hFigure);

% --- Check for Simulation Results ---
if isempty(handles.risultati_sim) || ~isstruct(handles.risultati_sim) || ~isfield(handles.risultati_sim, 'timeS') || isempty(handles.risultati_sim.timeS) % Check timeS
    disp('[Helper] plotResults: Cannot plot: no valid simulation results available.');
    return
end
res = handles.risultati_sim; % Uses names from runSimulation output

% --- Check for Stress Results (Optional) ---
stressRes = []; % Initialize
if isfield(handles, 'stressResults') && isstruct(handles.stressResults) && isfield(handles.stressResults, 'timeS') && ~isempty(handles.stressResults.timeS) && length(handles.stressResults.timeS) == length(res.timeS)
    stressRes = handles.stressResults;
    disp('[Helper] plotResults: Found stress results to plot.');
else
    disp('[Helper] plotResults: Stress results not found or invalid. Skipping stress plots.');
end

% --- Check for Parameters ---
if isempty(handles.parametri) || ~isstruct(handles.parametri)
    disp('[Helper] plotResults: Cannot plot: parameters struct missing.');
    return;
end
params = handles.parametri; % Uses names from loadCaseParameters output

% --- Verify required fields (basic check) ---
% (Keep existing checks for res and params)
required_res = {'timeS', 'gasPressurePa', 'projectilePositionM', 'projectileVelocityMps', ...
                'angularVelocityRadps', 'remainingPropellantMassKg', 'gasMassKg', ...
                'frictionWorkJ', 'heatLossJ'};
 if ~all(isfield(res, required_res)) || length(res.timeS) < 2
     fprintf(2,'[Helper] plotResults: Error plotting: simulation results incomplete.\n');
     return;
 end
required_params = {'projMass_m', 'projMomentOfInertia_Ip', 'barrelLength_m', ...
                   'initialPropellantMass_m', 'impetus_F'};
 if ~all(isfield(params, required_params))
      missing_p = required_params(~isfield(params, required_params));
      fprintf(2,'[Helper] plotResults: Error plotting: parameters struct missing: %s\n', strjoin(missing_p,', '));
      return;
 end


disp('[Helper] plotResults: Plotting results...');
t_ms = res.timeS * 1000;
t_end_ms = t_ms(end);
endIndex = length(t_ms);

% --- Define Axes Handles and Layout ---
% Make space for 2 more plots (8 total) -> Use 4x2 layout
axesHandles = {handles.axesPressure, handles.axesVelocity, ... % Row 1
               handles.axesPosition, handles.axesOmega, ...    % Row 2
               handles.axesMass, handles.axesEnergy, ...       % Row 3
               [], []}; % Placeholders for Row 4 (Stress & Safety Factor)
axTags = {'axesPressure', 'axesVelocity', 'axesPosition', 'axesOmega', 'axesMass', 'axesEnergy', 'axesStress', 'axesSafetyFactor'};

% Check if Stress/SF axes exist, create if not (assuming rightPanel exists)
if ~isfield(handles,'axesStress') || ~ishandle(handles.axesStress)
    axMarginX = 0.05; axMarginY = 0.05; % Adjusted margins for 4 rows
    axWidth = (1 - 3*axMarginX) / 2;
    axHeight = (1 - 5*axMarginY) / 4; % Height for 4 rows
    xPos1 = axMarginX; xPos2 = 2*axMarginX + axWidth;
    yPos4 = axMarginY; % Bottom row Y position

    handles.axesStress = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos1, yPos4, axWidth, axHeight], 'Tag', 'axesStress');
    handles.axesSafetyFactor = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos2, yPos4, axWidth, axHeight], 'Tag', 'axesSafetyFactor');
    axesHandles{7} = handles.axesStress;
    axesHandles{8} = handles.axesSafetyFactor;
    guidata(hFigure, handles); % Save new handles
    disp('Created axes for Stress and Safety Factor plots.');
else
    axesHandles{7} = handles.axesStress;
    axesHandles{8} = handles.axesSafetyFactor;
    disp('Using existing axes for Stress and Safety Factor plots.');
end


try
    % --- Pressure Plot ---
    hAx = handles.axesPressure; if ~ishandle(hAx), error('Invalid axesPressure handle'); end
    cla(hAx); hold(hAx, 'on');
    % (Keep existing pressure plotting code: shaded zones, main line, max marker)
    ax_limits_p = [0, t_end_ms, 0, max(res.gasPressurePa / 1e6) * 1.1]; axis(hAx, ax_limits_p);
    x_min_p = ax_limits_p(1); x_max_p = ax_limits_p(2); y_max_p_est = ax_limits_p(4);
    yellow_low = 390; yellow_high = 415; red_low = 415;
    yellow_color = [1 1 0.6]; red_color = [1 0.4 0.4]; alpha_value = 0.3;
    x_vertices_p = [x_min_p, x_max_p, x_max_p, x_min_p];
    y_vertices_yellow = [yellow_low, yellow_low, yellow_high, yellow_high];
    if yellow_high > ax_limits_p(3) && yellow_low < y_max_p_est, patch(hAx, x_vertices_p, y_vertices_yellow, yellow_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value, 'HandleVisibility', 'off'); end
    red_zone_top = max(red_low, y_max_p_est); y_vertices_red = [red_low, red_low, red_zone_top, red_zone_top];
    if red_zone_top > ax_limits_p(3) && red_low < y_max_p_est && red_zone_top > red_low, patch(hAx, x_vertices_p, y_vertices_red, red_color, 'EdgeColor', 'none', 'FaceAlpha', alpha_value, 'HandleVisibility', 'off'); end
    plot(hAx, t_ms, res.gasPressurePa / 1e6, 'LineWidth', 1.5, 'DisplayName', 'Pressure');
    pressureMPa = res.gasPressurePa / 1e6; [maxP_MPa, idxP] = max(pressureMPa);
    if ~isempty(idxP), timeP_max_ms = t_ms(idxP(1)); plot(hAx, timeP_max_ms, maxP_MPa, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Max Pressure'); text(hAx, timeP_max_ms, maxP_MPa, sprintf(' %.1f MPa', maxP_MPa), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold'); end
    grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Pressure [MPa]'); title(hAx, 'Pressure'); xlim(hAx, [0, t_end_ms]); axis(hAx,'tight'); ylim(hAx, [0 max(ylim(hAx))*1.05]); hold(hAx, 'off'); legend(hAx, 'show', 'Location', 'best');

    % --- Velocity Plot ---
    hAx = handles.axesVelocity; if ~ishandle(hAx), error('Invalid axesVelocity handle'); end
    cla(hAx); hold(hAx, 'on');
    % (Keep existing velocity plotting code: main line, max marker)
    plot(hAx, t_ms, res.projectileVelocityMps, 'LineWidth', 1.5, 'DisplayName', 'Velocity');
    velocityMps = res.projectileVelocityMps; [maxV_mps, idxV] = max(velocityMps);
    if ~isempty(idxV), timeV_max_ms = t_ms(idxV(1)); plot(hAx, timeV_max_ms, maxV_mps, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Max Velocity'); text(hAx, timeV_max_ms, maxV_mps, sprintf('%.1f m/s ', maxV_mps), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold'); end
    grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Velocity [m/s]'); title(hAx, 'Velocity'); xlim(hAx, [0, t_end_ms]); ylim(hAx, [0 max(ylim(hAx))*1.05]); hold(hAx, 'off'); legend(hAx, 'Location', 'best');

    % --- Position Plot ---
    hAx = handles.axesPosition; if ~ishandle(hAx), error('Invalid axesPosition handle'); end
    cla(hAx); hold(hAx, 'on');
    % (Keep existing position plotting code: trajectory, end marker, barrel length line)
    plot(hAx, t_ms, res.projectilePositionM * 100, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');
    plot(hAx, t_ms(endIndex), res.projectilePositionM(endIndex) * 100, 'r*', 'MarkerSize', 8, 'DisplayName', 'End Sim.');
    barrelLenParam = params.barrelLength_m;
    if ~isnan(barrelLenParam) && isnumeric(barrelLenParam), line(hAx, xlim(hAx), [barrelLenParam * 100, barrelLenParam * 100], 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'Barrel Length'); end
    hold(hAx, 'off'); grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Position [cm]'); title(hAx, 'Projectile Position'); legend(hAx, 'Location','northwest'); xlim(hAx, [0, t_end_ms]); ylim(hAx, [0 max(ylim(hAx))*1.05]);

    % --- Omega Plot ---
    hAx = handles.axesOmega; if ~ishandle(hAx), error('Invalid axesOmega handle'); end
    cla(hAx); hold(hAx, 'on');
    % (Keep existing omega plotting code: main line, max marker)
    omegaRpm = res.angularVelocityRadps * (60 / (2 * pi));
    plot(hAx, t_ms, omegaRpm, 'LineWidth', 1.5, 'DisplayName', 'Angular Velocity');
    [maxOmega_rpm, idxOmega] = max(omegaRpm);
    if ~isempty(idxOmega), timeOmega_max_ms = t_ms(idxOmega(1)); plot(hAx, timeOmega_max_ms, maxOmega_rpm, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Max RPM'); text(hAx, timeOmega_max_ms, maxOmega_rpm, sprintf('%.0f RPM ', maxOmega_rpm), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold'); end
    grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Velocity [RPM]'); title(hAx, 'Angular Velocity'); xlim(hAx, [0, t_end_ms]); ylim(hAx, [0 max(ylim(hAx))*1.05]); hold(hAx, 'off'); legend(hAx, 'Location', 'best');

    % --- Mass Plot ---
    hAx = handles.axesMass; if ~ishandle(hAx), error('Invalid axesMass handle'); end
    cla(hAx); hold(hAx, 'on');
    % (Keep existing mass plotting code: solid, gas, burnout line)
    plot(hAx, t_ms, res.remainingPropellantMassKg * 1000, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Solid');
    plot(hAx, t_ms, res.gasMassKg * 1000, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Gas');
    initialMassKg = params.initialPropellantMass_m;
    if initialMassKg > 1e-9, burnoutThresholdKg = 0.05 * initialMassKg; idxBurnout = find(res.remainingPropellantMassKg <= burnoutThresholdKg, 1, 'first'); if ~isempty(idxBurnout), timeBurnout_ms = t_ms(idxBurnout); yLimsMass = ylim(hAx); hLineBurnout = line(hAx, [timeBurnout_ms timeBurnout_ms], yLimsMass, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', '95% Burnout'); uistack(hLineBurnout,'bottom'); end; end
    hold(hAx, 'off'); grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Mass [g]'); title(hAx, 'Masses'); legend(hAx, 'Location', 'best'); xlim(hAx, [0, t_end_ms]); ylim(hAx, [0 max(ylim(hAx))*1.05]);

    % --- Energy/Work Plot ---
    hAx = handles.axesEnergy; if ~ishandle(hAx), error('Invalid axesEnergy handle'); end
    % (Keep existing energy plotting code: KE lin, KE rot, Work Bore Res, Q Loss)
    KE_lin = 0.5 * params.projMass_m * res.projectileVelocityMps.^2;
    KE_rot = 0.5 * params.projMomentOfInertia_Ip * res.angularVelocityRadps.^2;
    cla(hAx); hold(hAx, 'on');
    plot(hAx, t_ms, KE_lin / 1000, 'Color', [0 0.6 0], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'KE Linear');
    plot(hAx, t_ms, KE_rot / 1000, 'Color', [0.8 0.6 0], 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'KE Rotational');
    plot(hAx, t_ms, res.frictionWorkJ / 1000, 'Color', [0.8 0 0], 'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', 'Work Bore Res.');
    plot(hAx, t_ms, res.heatLossJ / 1000, 'Color', [0 0 0.8], 'LineStyle', ':', 'LineWidth', 1.5, 'DisplayName', 'Q Loss');
    hold(hAx, 'off'); grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Energy [kJ]'); title(hAx, 'Energy / Work'); legend(hAx, 'Location', 'best'); xlim(hAx, [0, t_end_ms]); ylim(hAx, [0 max(ylim(hAx))*1.05]);

    % --- NEW: Stress Plot ---
    hAx = handles.axesStress; if ~ishandle(hAx), error('Invalid axesStress handle'); end
    cla(hAx); % Clear axes
    if ~isempty(stressRes)
        hold(hAx, 'on');
        plot(hAx, t_ms, stressRes.sigma_theta_inner_Pa / 1e6, 'b-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{\theta} (Hoop)');
        plot(hAx, t_ms, stressRes.sigma_r_inner_Pa / 1e6, 'r--', 'LineWidth', 1.5, 'DisplayName', '\sigma_{r} (Radial)');
        plot(hAx, t_ms, stressRes.sigma_eq_GT_inner_Pa / 1e6, 'k-.', 'LineWidth', 1.5, 'DisplayName', '\sigma_{eq} (Tresca)');
        % Plot Yield Strength Line
        yline(hAx, stressRes.yield_strength_Pa / 1e6, 'm:', 'LineWidth', 1, 'DisplayName', 'Yield Strength');
        hold(hAx, 'off');
        grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Stress [MPa]');
        title(hAx, 'Stresses at Inner Barrel Wall'); legend(hAx, 'show', 'Location', 'best');
        xlim(hAx, [0, t_end_ms]); ylim(hAx, 'auto'); % Let y-axis scale automatically
    else
        title(hAx, 'Stress Results (N/A)'); axis(hAx, 'off');
    end

    % --- NEW: Safety Factor Plot ---
    hAx = handles.axesSafetyFactor; if ~ishandle(hAx), error('Invalid axesSafetyFactor handle'); end
    cla(hAx); % Clear axes
    if ~isempty(stressRes)
        hold(hAx, 'on');
        plot(hAx, t_ms, stressRes.safety_factor_inner, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Safety Factor');
        % Plot Minimum Safety Factor line and marker
        yline(hAx, stressRes.min_safety_factor, 'r--', 'LineWidth', 1, 'DisplayName', sprintf('Min SF = %.2f', stressRes.min_safety_factor));
        plot(hAx, stressRes.time_at_min_sf_s * 1000, stressRes.min_safety_factor, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'HandleVisibility', 'off');
        yline(hAx, 1.0, 'k:', 'LineWidth', 1, 'DisplayName', 'SF = 1.0'); % Indicate failure threshold
        hold(hAx, 'off');
        grid(hAx, 'on'); xlabel(hAx, 'Time [ms]'); ylabel(hAx, 'Safety Factor [-]');
        title(hAx, 'Safety Factor at Inner Wall'); legend(hAx, 'show', 'Location', 'best');
        xlim(hAx, [0, t_end_ms]);
        % Set Y limits appropriately, e.g., start near 0 or min_sf, maybe cap max SF
        min_plot_sf = max(0, stressRes.min_safety_factor - 0.5);
        max_plot_sf = max(stressRes.min_safety_factor*1.5, 2.0); % Adjust upper limit as needed
        ylim(hAx, [min_plot_sf, max_plot_sf]);
    else
        title(hAx, 'Safety Factor Results (N/A)'); axis(hAx, 'off');
    end


    datacursormode(hFigure, 'on');
    disp('[Helper] plotResults: Plotting completed (including stress plots if available).');

catch ME
     fprintf(2,'[Helper] plotResults: ERROR during plotting: %s\n', ME.message);
     if ~isempty(ME.stack)
        fprintf(2,'File: %s, Line: %d\n', ME.stack(1).file, ME.stack(1).line);
     end
     try datacursormode(hFigure, 'off'); catch; end
end

end % End function plotResults