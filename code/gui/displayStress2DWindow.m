% --- gui/displayStress2DWindow.m ---
% MODIFIED: Added colored 2D cross-section stress plot.
% CORRECTED: Fixed missing 'end' for inner try-catch block.
% CORRECTED: Renamed function to match filename convention.
% REMOVED: Custom blue-red colormap definition (using 'jet').
function displayStress2DWindow(stressResults, barrelData)
% Creates a new window showing 2D stress & displacement distribution plots
% and a colored 2D cross-section view of equivalent stress.

    windowTitle = 'Barrel Stress & Displacement Analysis'; 
    
    % --- Input Validation ---
    if nargin < 2 || isempty(stressResults) || isempty(barrelData)
        warning('%s: Missing or empty input data.', windowTitle);
        return;
    end
    
    % --- Check if necessary radial data exists ---
    required_radial_fields = {'radius_vector_m', 'sigma_r_radial_Pa', ...
                              'sigma_theta_radial_Pa', 'sigma_eq_GT_radial_Pa', ...
                              'pressure_at_max_stress_Pa', 'yield_strength_Pa', ...
                              'radial_displacement_m'}; 
    if ~all(isfield(stressResults, required_radial_fields))
         warning('%s: stressResults struct missing required radial data fields.', windowTitle);
         return;
    end
    
    % --- Check material properties for displacement ---
    can_calc_displacement = false;
    if isfield(barrelData, 'youngsModulus_Pa') && isfield(barrelData, 'poissonsRatio')
        E_mod = barrelData.youngsModulus_Pa;
        nu = barrelData.poissonsRatio;
        if ~isnan(E_mod) && isnumeric(E_mod) && E_mod > 0 && ~isnan(nu) && isnumeric(nu) && nu >= 0 && nu < 0.5
            can_calc_displacement = true;
            stressResults.youngsModulus_Pa = E_mod; % Store valid props used
            stressResults.poissonsRatio = nu;
        else
             warning('displayStress2DWindow: Invalid Young''s Modulus (%.2e Pa) or Poisson''s Ratio (%.2f) in barrelData. Displacement calc skipped.', E_mod, nu);
        end
    else
         warning('displayStress2DWindow: Missing Young''s Modulus or Poisson''s Ratio in barrelData. Displacement calc skipped.');
    end

    % --- Extract data ---
    radius_m = stressResults.radius_vector_m;
    sigma_r_MPa = stressResults.sigma_r_radial_Pa / 1e6;
    sigma_t_MPa = stressResults.sigma_theta_radial_Pa / 1e6;
    sigma_eq_MPa = stressResults.sigma_eq_GT_radial_Pa / 1e6; 
    yield_MPa = stressResults.yield_strength_Pa / 1e6;
    p_max_stress_MPa = stressResults.pressure_at_max_stress_Pa / 1e6;
    u_r_m = stressResults.radial_displacement_m; 
    r_i = barrelData.boreDiameter_m / 2; 
    r_e = barrelData.outerDiameter_m / 2; 
    
    % --- Create Figure ---
    hFig = figure('Name', windowTitle, ...
                  'Units', 'normalized', 'Position', [0.15 0.1 0.7 0.8], ... 
                  'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure');
                  
    % --- Plot 1: Stress Distribution vs Radius ---
    hAx1 = subplot(3, 1, 1, 'Parent', hFig); 
    hold(hAx1, 'on');
    plot(hAx1, radius_m * 1000, sigma_t_MPa, 'b-', 'LineWidth', 1.5, 'DisplayName', '\sigma_{\theta} (Hoop)');
    plot(hAx1, radius_m * 1000, sigma_r_MPa, 'r--', 'LineWidth', 1.5, 'DisplayName', '\sigma_{r} (Radial)');
    plot(hAx1, radius_m * 1000, sigma_eq_MPa, 'k-.', 'LineWidth', 1.5, 'DisplayName', '\sigma_{eq} (Tresca)');
    yline(hAx1, yield_MPa, 'm:', 'LineWidth', 1, 'DisplayName', 'Yield Strength');
    hold(hAx1, 'off');
    grid(hAx1, 'on'); xlabel(hAx1, 'Radius [mm]'); ylabel(hAx1, 'Stress [MPa]');
    titleStr1 = sprintf('Stress Distribution @ Max Eq. Stress (P = %.1f MPa)', p_max_stress_MPa);
    title(hAx1, titleStr1); legend(hAx1, 'show', 'Location', 'best');
    xlim(hAx1, [r_i*1000, r_e*1000]); ylim(hAx1, 'auto');
    
    % --- Plot 2: Radial Displacement vs Radius ---
    hAx2 = subplot(3, 1, 2, 'Parent', hFig); 
    if can_calc_displacement && ~all(isnan(u_r_m)) % Check again if displacement is valid
        plot(hAx2, radius_m * 1000, u_r_m * 1e6, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Radial Displacement'); % Plot in micrometers (um)
        grid(hAx2, 'on'); xlabel(hAx2, 'Radius [mm]'); ylabel(hAx2, 'Displacement [\mum]'); 
        % Use stored E and nu if they were valid
        E_GPa_disp = stressResults.youngsModulus_Pa / 1e9;
        nu_disp = stressResults.poissonsRatio;
        titleStr2 = sprintf('Radial Displacement @ Max Eq. Stress (E=%.0f GPa, \\nu=%.2f)', E_GPa_disp, nu_disp);
        title(hAx2, titleStr2); legend(hAx2, 'show', 'Location', 'best');
        xlim(hAx2, [r_i*1000, r_e*1000]); ylim(hAx2, 'auto');
    else
        title(hAx2, 'Radial Displacement (N/A - Check Material Properties)');
        axis(hAx2, 'off');
    end
    
    % --- Plot 3: Colored 2D Cross-Section ---
    hAx3 = subplot(3, 1, 3, 'Parent', hFig); 
    try % OUTER TRY
        % Create a grid for the cross-section plot
        n_grid = 200; 
        max_rad_plot = r_e * 1.05; 
        x_lin = linspace(-max_rad_plot, max_rad_plot, n_grid);
        y_lin = linspace(-max_rad_plot, max_rad_plot, n_grid);
        [X, Y] = meshgrid(x_lin, y_lin);
        R = sqrt(X.^2 + Y.^2); 
        
        % Interpolate the equivalent stress onto the grid
        [unique_radius_m, ia, ~] = unique(radius_m);
        unique_sigma_eq_MPa = sigma_eq_MPa(ia);
        Stress_interp = interp1(unique_radius_m, unique_sigma_eq_MPa, R, 'linear', NaN); 
        
        % Mask out areas outside the barrel material
        Stress_interp(R < r_i | R > r_e) = NaN;
        
        % Plot using pcolor
        pcolor(hAx3, X * 1000, Y * 1000, Stress_interp); 
        shading(hAx3, 'interp'); 
        axis(hAx3, 'equal'); 
        colormap(hAx3, 'jet'); % Use jet colormap
        hcb = colorbar(hAx3);
        ylabel(hcb, 'Equivalent Stress (Tresca) [MPa]');
        xlabel(hAx3, 'X [mm]'); ylabel(hAx3, 'Y [mm]');
        title(hAx3, 'Equivalent Stress Distribution in Cross-Section (@ Max)');
        
        % Set color limits 
        % Calculate max stress, handling potential all-NaN case
        max_stress_plot = max(Stress_interp(:));
        if isempty(max_stress_plot) || isnan(max_stress_plot)
             max_stress_plot = yield_MPa; % Fallback if no valid stress data
        end
        clim(hAx3, [0, max(max_stress_plot, yield_MPa)]); % Use clim (newer)
        
        grid(hAx3, 'off'); 

        % --- Aggiungi testo con valore massimo stress ---
        try % INNER TRY for text annotation
            max_stress_value_MPa = max(Stress_interp(~isnan(Stress_interp)));             
            if ~isempty(max_stress_value_MPa) && isfinite(max_stress_value_MPa)
                max_stress_text = sprintf('Max Eq. Stress: %.1f MPa', max_stress_value_MPa);
                text(hAx3, 0.05, 0.95, max_stress_text, ...
                     'Units', 'normalized', ... 
                     'VerticalAlignment', 'top', ...
                     'HorizontalAlignment', 'left', ...
                     'Color', 'white', ... 
                     'BackgroundColor', [0 0 0 0.5], ... 
                     'FontSize', 9, ...
                     'FontWeight', 'bold');
            end
        catch ME_text % INNER CATCH
            warning('Could not add max stress text annotation: %s', ME_text.message);
        end % <<<=== 'end' CORRETTO PER IL BLOCCO try...catch ME_text ===>>>
        % --- Fine aggiunta testo ---

    catch ME_plot_section % OUTER CATCH 
        title(hAx3, 'Error plotting 2D section');
        axis(hAx3,'on'); grid(hAx3,'on');
        warning('Error creating 2D stress section plot: %s', ME_plot_section.message);
    end % OUTER END

end % --- End displayStress2DWindow ---
