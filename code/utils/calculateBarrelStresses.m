% utils/calculateBarrelStresses.m
% MODIFIED: Calculates stress AND radial displacement distribution
%           across radius at time of max eq. stress.
function stressResults = calculateBarrelStresses(simulationResults, barrelData)
% Calculates stresses/displacement at inner wall AND distributions
% across radius at the time of maximum equivalent stress.
% Uses Lamé's equations & Guest-Tresca.
%
% INPUTS:
%   simulationResults: Struct with sim results (timeS, gasPressurePa).
%   barrelData: Struct with barrel params (boreDiameter_m, outerDiameter_m,
%               yieldStrength_Pa, youngsModulus_Pa, poissonsRatio). % Added E, nu
%
% OUTPUT:
%   stressResults: Struct containing time series and radial distribution results:
%                  - timeS, pressurePa, sigma_r_inner_Pa, sigma_theta_inner_Pa,
%                  - sigma_eq_GT_inner_Pa, safety_factor_inner
%                  - min_safety_factor, time_at_min_sf_s
%                  - max_eq_stress_Pa, time_at_max_stress_s
%                  - yield_strength_Pa
%                  - pressure_at_max_stress_Pa
%                  - radius_vector_m
%                  - sigma_r_radial_Pa, sigma_theta_radial_Pa, sigma_eq_GT_radial_Pa
%                  - youngsModulus_Pa, poissonsRatio % Store used material props
%                  - radial_displacement_m: Radial displacement vs radius at max stress time [m] % NEW

    disp('Calculating barrel stresses & displacement (inner wall and radial @ max)...');
    stressResults = struct(); % Initialize output

    % --- Input Validation ---
    if ~isstruct(simulationResults) || ~isfield(simulationResults, 'timeS') || ~isfield(simulationResults, 'gasPressurePa') || isempty(simulationResults.timeS)
        error('calculateBarrelStresses: Invalid simulationResults input.');
    end
    % Check for material properties needed for displacement
    required_barrel_fields = {'boreDiameter_m', 'outerDiameter_m', 'yieldStrength_Pa', 'youngsModulus_Pa', 'poissonsRatio'};
    missing_fields = required_barrel_fields(~isfield(barrelData, required_barrel_fields));
    if ~isempty(missing_fields)
         warning('calculateBarrelStresses: barrelData missing fields required for stress/displacement (%s). Displacement calc skipped.', strjoin(missing_fields,', '));
         can_calc_displacement = false;
    else
         % Check if values are valid
         E_mod = barrelData.youngsModulus_Pa;
         nu = barrelData.poissonsRatio;
         if isnan(E_mod) || ~isnumeric(E_mod) || E_mod <= 0 || isnan(nu) || ~isnumeric(nu) || nu < 0 || nu >= 0.5
             warning('calculateBarrelStresses: Invalid Young''s Modulus (%.2e Pa) or Poisson''s Ratio (%.2f) in barrelData. Displacement calc skipped.', E_mod, nu);
             can_calc_displacement = false;
         else
             can_calc_displacement = true;
             stressResults.youngsModulus_Pa = E_mod; % Store valid props used
             stressResults.poissonsRatio = nu;
         end
    end


    % Extract data
    t = simulationResults.timeS(:);
    p_i = simulationResults.gasPressurePa(:); p_i(p_i < 0) = 0;
    r_i = barrelData.boreDiameter_m / 2; r_e = barrelData.outerDiameter_m / 2;
    sigma_y = barrelData.yieldStrength_Pa;

    if r_i <= 0 || r_e <= r_i || sigma_y <= 0
        error('calculateBarrelStresses: Invalid barrel dimensions or yield strength.');
    end

    % --- Calculate Stresses at Inner Wall (Time Series) ---
    a = r_e / r_i; a2 = a^2;
    if abs(a2 - 1) < 1e-9, error('calculateBarrelStresses: Inner and outer radii too close.'); end

    sigma_r_inner_Pa = -p_i;
    sigma_theta_inner_Pa = p_i .* (a2 + 1) / (a2 - 1);
    sigma_eq_GT_inner_Pa = p_i .* (2 * a2) / (a2 - 1);

    safety_factor_inner = zeros(size(p_i));
    nonZeroStressIdx = sigma_eq_GT_inner_Pa > 1e-3;
    safety_factor_inner(nonZeroStressIdx) = sigma_y ./ sigma_eq_GT_inner_Pa(nonZeroStressIdx);
    safety_factor_inner(~nonZeroStressIdx) = Inf;

    % --- Find Max Stress and Min Safety Factor (Time Series) ---
    [min_sf, min_sf_idx] = min(safety_factor_inner); time_at_min_sf = t(min_sf_idx);
    [max_stress, max_stress_idx] = max(sigma_eq_GT_inner_Pa); time_at_max_stress = t(max_stress_idx);
    pressure_at_max_stress = p_i(max_stress_idx);

    % --- Store Time Series Results ---
    stressResults.timeS = t; stressResults.pressurePa = p_i;
    stressResults.sigma_r_inner_Pa = sigma_r_inner_Pa; stressResults.sigma_theta_inner_Pa = sigma_theta_inner_Pa;
    stressResults.sigma_eq_GT_inner_Pa = sigma_eq_GT_inner_Pa; stressResults.safety_factor_inner = safety_factor_inner;
    stressResults.min_safety_factor = min_sf; stressResults.time_at_min_sf_s = time_at_min_sf;
    stressResults.max_eq_stress_Pa = max_stress; stressResults.time_at_max_stress_s = time_at_max_stress;
    stressResults.yield_strength_Pa = sigma_y; stressResults.pressure_at_max_stress_Pa = pressure_at_max_stress;

    % --- Calculate Radial Stress & Displacement Distribution at Max Stress Time ---
    fprintf('Calculating radial stress/displacement distribution at P = %.2f MPa (t=%.4f ms)...\n', pressure_at_max_stress/1e6, time_at_max_stress*1000);
    num_radial_points = 100;
    r_vector = linspace(r_i, r_e, num_radial_points); % Radius vector [m]
    r2_vector = r_vector.^2;
    p_load = pressure_at_max_stress;
    lame_const_factor = (p_load * r_i^2) / (r_e^2 - r_i^2);

    % Stresses
    sigma_r_radial = lame_const_factor * (1 - (r_e^2 ./ r2_vector));
    sigma_theta_radial = lame_const_factor * (1 + (r_e^2 ./ r2_vector));
    sigma_eq_GT_radial = sigma_theta_radial - sigma_r_radial;

    % Displacement (u) - Initialize to NaN
    radial_displacement = nan(size(r_vector));
    if can_calc_displacement
        E_mod = stressResults.youngsModulus_Pa;
        nu = stressResults.poissonsRatio;
        disp_const_factor = (p_load * r_i^2) / (E_mod * (r_e^2 - r_i^2));
        % Formula for radial displacement u(r)
        radial_displacement = disp_const_factor * ((1 - nu) * r_vector + (1 + nu) * (r_e^2 ./ r_vector));
    end

    % Store radial results
    stressResults.radius_vector_m = r_vector;
    stressResults.sigma_r_radial_Pa = sigma_r_radial;
    stressResults.sigma_theta_radial_Pa = sigma_theta_radial;
    stressResults.sigma_eq_GT_radial_Pa = sigma_eq_GT_radial;
    stressResults.radial_displacement_m = radial_displacement; % Store displacement (may be NaN)
    % --- END RADIAL CALCULATION SECTION ---

    fprintf('Barrel stress/displacement calculation complete. Min Safety Factor = %.2f\n', min_sf);

end