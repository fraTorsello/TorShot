% gui/plotBarrel3D.m
function plotBarrel3D(hAxes, barrelData, stressValuePa, yieldStrengthPa)
% Plots a simplified 3D model of the barrel on the specified axes.
% Optionally colors the outer surface based on a provided stress value
% relative to the yield strength.
% Assumes barrelData contains outerDiameter_m, boreDiameter_m, barrelLength_m.
% INPUTS:
%   hAxes:           Handle to the axes where the plot should be drawn.
%   barrelData:      Struct containing barrel parameters in SI units.
%   stressValuePa:   (Optional) Stress value [Pa] to use for coloring the
%                    outer surface (e.g., equivalent stress at max pressure).
%   yieldStrengthPa: (Optional) Yield strength [Pa] used as the upper limit
%                    for the color scale. Required if stressValuePa is provided.

    plotWithStress = false;
    if nargin >= 4 && isnumeric(stressValuePa) && ~isnan(stressValuePa) && isnumeric(yieldStrengthPa) && ~isnan(yieldStrengthPa) && yieldStrengthPa > 0
        plotWithStress = true;
    elseif nargin >= 3 && ( ~isnumeric(stressValuePa) || isnan(stressValuePa) )
        % If stressValuePa provided but invalid, ignore it
         fprintf('plotBarrel3D: Invalid stressValuePa provided. Plotting without stress colors.\n');
    elseif nargin == 3 % Stress provided but yield strength missing
         warning('plotBarrel3D: stressValuePa provided but yieldStrengthPa missing or invalid. Plotting without stress colors.');
    end % else nargin < 3, plot basic model


    if nargin < 2 || ~ishandle(hAxes) || ~isstruct(barrelData)
        warning('plotBarrel3D: Invalid input axes handle or barrelData struct.');
        if ishandle(hAxes), cla(hAxes); title(hAxes, 'Error: Invalid Data'); axis(hAxes,'off'); end
        return;
    end

    % --- Extract required parameters ---
    try
        outerDia = barrelData.outerDiameter_m; boreDia = barrelData.boreDiameter_m; len = barrelData.barrelLength_m; barrelName = barrelData.barrelName;
        if isnan(outerDia) || isnan(boreDia) || isnan(len) || outerDia <= 0 || boreDia <= 0 || len <= 0 || boreDia >= outerDia, error('Invalid dimensions'); end
    catch ME_extract, warning('plotBarrel3D: Error extracting parameters: %s', ME_extract.message); cla(hAxes); title(hAxes, 'Error: Missing/Invalid Dimensions'); axis(hAxes,'off'); return; end

    % --- Plotting Setup ---
    n_surf = 50; r_outer = outerDia / 2; r_bore = boreDia / 2;
    [Xo, Yo, Zo] = cylinder(r_outer, n_surf); Zo = Zo * len;
    [Xi, Yi, Zi] = cylinder(r_bore, n_surf); Zi = Zi * len;

    % --- Clear and Plot on specified axes ---
    cla(hAxes); hold(hAxes, 'on');

    % --- Plot Surfaces ---
    hSurfOuter = []; hSurfBore = []; % Initialize handles
    if plotWithStress
        % Create color data based on the single stress value
        % Option 1: Uniform color based on the single value
        % colorLevel = min(1, max(0, stressValuePa / yieldStrengthPa)); % Normalize 0-1
        % C_outer = ones(size(Zo)) * colorLevel; % Apply uniform color level

        % Option 2: Create a gradient (more complex, let's stick to uniform for now)
        % For uniform color, directly set FaceColor based on colormap lookup
        cmap = jet(256); % Use jet colormap
        stressRatio = min(1, max(0, stressValuePa / yieldStrengthPa));
        colorIndex = round(1 + stressRatio * (size(cmap, 1) - 1));
        faceColor = cmap(colorIndex, :);
        hSurfOuter = surf(hAxes, Xo, Yo, Zo, 'FaceColor', faceColor, 'EdgeColor', 'none', 'DisplayName', 'Outer Surface (Stress)');

        % Add colorbar
        colormap(hAxes, cmap);
        hcb = colorbar(hAxes);
        ylabel(hcb, 'Equiv. Stress / Yield Strength');
        caxis(hAxes, [0, 1]); % Set color limits from 0 to 1 (normalized)
        % Alternative: caxis(hAxes, [0, yieldStrengthPa/1e6]); ylabel(hcb, 'Eq. Stress (MPa)');

    else % Plot without stress colors
        hSurfOuter = surf(hAxes, Xo, Yo, Zo, 'FaceColor', [0.6 0.6 0.65], 'EdgeColor', 'none', 'DisplayName', 'Outer Surface');
        % Remove colorbar if it exists from a previous plot
        hcb_old = findobj(gcf, 'Type', 'colorbar', 'Parent', hAxes); % Find colorbar associated with hAxes
        delete(hcb_old);
        colormap(hAxes, gray); % Reset colormap
    end

    hSurfBore = surf(hAxes, Xi, Yi, Zi, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'DisplayName', 'Bore');

    % Plot End Caps (same as before)
    theta_cap = linspace(0, 2*pi, n_surf+1); x_cap_outer = r_outer * cos(theta_cap); y_cap_outer = r_outer * sin(theta_cap);
    x_cap_inner = r_bore * cos(theta_cap); y_cap_inner = r_bore * sin(theta_cap);
    patch(hAxes, x_cap_outer, y_cap_outer, zeros(size(x_cap_outer)), [0.6 0.6 0.65], 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Breech Face', 'FaceAlpha', 0.8);
    patch(hAxes, x_cap_inner, y_cap_inner, zeros(size(x_cap_inner)), [1 1 1], 'EdgeColor', 'k', 'LineWidth', 0.5, 'FaceAlpha', 0.8);
    patch(hAxes, x_cap_outer, y_cap_outer, ones(size(x_cap_outer))*len, [0.6 0.6 0.65], 'EdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Muzzle Face', 'FaceAlpha', 0.8);
    patch(hAxes, x_cap_inner, y_cap_inner, ones(size(x_cap_inner))*len, [1 1 1], 'EdgeColor', 'k', 'LineWidth', 0.5, 'FaceAlpha', 0.8);

    hold(hAxes, 'off');

    % --- Set Axes Properties ---
    axis(hAxes, 'equal'); xlabel(hAxes,'X (m)'); ylabel(hAxes,'Y (m)'); zlabel(hAxes,'Z (Length, m)');
    titleStr = {barrelName, '(Simplified Model)'};
    if plotWithStress, titleStr{2} = [titleStr{2}, ' - Stress Color @ Max P']; end
    title(hAxes, titleStr, 'Interpreter', 'none');
    grid(hAxes, 'on'); view(hAxes, 3); camlight(hAxes, 'headlight'); material(hAxes, 'dull');
    max_radius = outerDia/2; axis(hAxes, [-max_radius*1.1, max_radius*1.1, -max_radius*1.1, max_radius*1.1, -len*0.05, len*1.05]);

end