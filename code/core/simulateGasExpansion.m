% --- MODIFIED simulateGasExpansion.m ---
% Calculates h_conv dynamically based on pressure

function [expansionResults, finalGasEnergy] = simulateGasExpansion(initialState, constantGasMass, parameters, barrelParameters) % REMOVED expansionHeatCoeff_h argument
%simulateGasExpansion Simulates the gas expansion phase after propellant burnout.
%   MODIFIED: Calculates h_conv dynamically based on pressure.
%   MODIFIED: Simulation continues until maxSafetyTime_s or until
%             gas pressure drops to ambientPressure.
%             Removed the barrel exit event specific to this phase.
%
%   [expansionResults, finalGasEnergy] = simulateGasExpansion(initialState, constantGasMass, parameters, barrelParameters) % UPDATED Signature
%   Simulates gas expansion and projectile motion starting from a defined
%   initial state, considering expansion work and heat transfer to the
%   barrel (using a dynamically calculated heat transfer coefficient).
%   Gas mass is assumed constant during this phase.
%   Friction is neglected in this simplified version.
%
%   INPUTS:
%   initialState: Struct with the initial system state. Fields: t_start, T_gas0, x_p0, v_p0
%   constantGasMass: Total gas mass [kg]
%   parameters: Struct with general parameters. Required fields:
%               specificGasConstant_R, specificHeatRatio_gamma, covolume_b,
%               projMass_m, projArea_Ab, boreDiameter_Db, initialTemperature_K,
%               maxSafetyTime_s, barrelLength_m, initialFreeVolume_V0, ambientPressure,
%               h_conv_base_coefficient, h_conv_pressure_reference_Pa, h_conv_pressure_exponent % ADDED h_conv params
%   barrelParameters: Struct with barrel parameters (massKg, specificHeatJkgK) - (Only used for context now)
%   % REMOVED: expansionHeatCoeff_h argument
%
%   OUTPUTS:
%   expansionResults: Struct with expansion simulation results: t, T_gas, x_p, v_p, Q_loss_exp, P_gas_calc.
%   finalGasEnergy: Struct with final gas energies: U_gas_residual, KE_gas_estimated, Q_loss_during_expansion.

% --- Input Validation (Basic) ---
if ~isstruct(initialState) || ~isfield(initialState, 't_start') || ...
   ~isfield(initialState, 'T_gas0') || ~isfield(initialState, 'x_p0') || ...
   ~isfield(initialState, 'v_p0')
    error('The ''initialState'' structure is incomplete.');
end
if ~isnumeric(constantGasMass) || constantGasMass <= 0
    error('Gas mass must be a positive number.');
end

% MODIFIED: Check for base h_conv parameters instead of fixed h_conv
req_params = {'specificGasConstant_R', 'specificHeatRatio_gamma', 'covolume_b', ...
              'projMass_m', 'projArea_Ab', 'boreDiameter_Db', 'initialTemperature_K', ...
              'maxSafetyTime_s', 'barrelLength_m', 'initialFreeVolume_V0', 'ambientPressure', ...
              'h_conv_base_coefficient', 'h_conv_pressure_reference_Pa', 'h_conv_pressure_exponent'}; % Added h_conv params
missing_params = req_params(~isfield(parameters, req_params));
if ~isempty(missing_params)
     error('The ''parameters'' structure is missing required fields for expansion simulation: %s', strjoin(missing_params, ', '));
end
% Basic validation for new h_conv params
if ~isnumeric(parameters.h_conv_base_coefficient) || parameters.h_conv_base_coefficient < 0 || ...
   ~isnumeric(parameters.h_conv_pressure_reference_Pa) || parameters.h_conv_pressure_reference_Pa <= 0 || ...
   ~isnumeric(parameters.h_conv_pressure_exponent)
     error('Invalid h_conv base parameters found in the ''parameters'' struct.');
end

% --- Local Parameters ---
gasMass = constantGasMass;
specificGasConstant_R = parameters.specificGasConstant_R;
gamma = parameters.specificHeatRatio_gamma;
covolume_b = parameters.covolume_b;
projMass_m = parameters.projMass_m;
projArea_Ab = parameters.projArea_Ab;
boreDiameter_Db = parameters.boreDiameter_Db;
initialFreeVolume_V0 = parameters.initialFreeVolume_V0;
maxSafetyTime_s = parameters.maxSafetyTime_s;
ambientPressure = parameters.ambientPressure;

Cv_gas = specificGasConstant_R / (gamma - 1);
T_barrel = parameters.initialTemperature_K; % Barrel temperature assumed constant

% Extract h_conv base parameters
h_conv_base = parameters.h_conv_base_coefficient;
h_conv_pref = parameters.h_conv_pressure_reference_Pa;
h_conv_n    = parameters.h_conv_pressure_exponent;

% --- Initial Conditions for Expansion ODE ---
y0_exp = [initialState.T_gas0; initialState.x_p0; initialState.v_p0; 0]; % [T_gas; x_p; v_p; Q_loss_exp]

% --- Time Interval ---
t_start = initialState.t_start;
t_end_max = maxSafetyTime_s;
if t_start >= t_end_max
    error('Initial time t_start (%.4f s) is greater than or equal to maxSafetyTime_s (%.4f s)', t_start, t_end_max);
end
tspan_exp = [t_start, t_end_max];

% --- ODE Solver Options ---
options_exp = odeset(...
    'RelTol', 1e-5, ...
    'AbsTol', [1e-2, 1e-7, 1e-4, 1e-1], ... % Tolerances for [T, x, v, Q]
    'NonNegative', [1, 3, 4], ... % T, v, Q >= 0
    'Events', @(t,y) atmosphericPressureEvent(t, y, gasMass, specificGasConstant_R, covolume_b, initialFreeVolume_V0, projArea_Ab, ambientPressure) ...
);

% --- Call ODE Solver ---
fprintf('Starting expansion simulation from t=%.4f s (until P_amb=%.1f kPa or t=%.4f s)...\n', t_start, ambientPressure/1000, t_end_max);
try
    % MODIFIED: Pass h_conv base parameters to odesExpansionOnly instead of fixed h_conv
    [t, y_sol, te, ye, ie] = ode45(@(t,y) odesExpansionOnly(t, y, gasMass, Cv_gas, projMass_m, projArea_Ab, boreDiameter_Db, initialFreeVolume_V0, specificGasConstant_R, covolume_b, T_barrel, h_conv_base, h_conv_pref, h_conv_n), ...
                                   tspan_exp, y0_exp, options_exp);
    fprintf('Expansion simulation finished at t=%.4f s.\n', t(end));
    if ~isempty(ie)
         fprintf('  -> Atmospheric Pressure Event reached at t=%.6f s.\n', te(end));
    else
         warning('Expansion simulation stopped at maximum time t_max_safety (%.4f s) BEFORE reaching atmospheric pressure.', t(end));
    end
catch ME_ode_exp
    error('Error during expansion simulation: %s', ME_ode_exp.message);
end

% --- Organize Results ---
expansionResults = struct();
expansionResults.t = t;
expansionResults.T_gas = y_sol(:, 1);
expansionResults.x_p = y_sol(:, 2);
expansionResults.v_p = y_sol(:, 3);
expansionResults.Q_loss_exp = y_sol(:, 4); % Heat lost *during* this phase

% Recalculate Pressure for output (using Noble-Able)
P_gas_calc = zeros(size(t));
for i = 1:length(t)
    V_inst_i = initialFreeVolume_V0 + projArea_Ab * expansionResults.x_p(i);
    eff_vol_i = V_inst_i - gasMass * covolume_b;
    if eff_vol_i > 1e-9 % Avoid division by zero/small number
        P_gas_calc(i) = (gasMass * specificGasConstant_R * expansionResults.T_gas(i)) / eff_vol_i;
    else
        P_gas_calc(i) = P_gas_calc(max(1, i-1)); % Maintain previous pressure if volume invalid
    end
end
P_gas_calc = max(ambientPressure, P_gas_calc); % Ensure pressure >= ambient
expansionResults.P_gas_calc = P_gas_calc;


% --- Calculate Final Gas Energies ---
finalGasEnergy = struct();
T_gas_final = expansionResults.T_gas(end);
v_proj_final = expansionResults.v_p(end); % Velocity at the final instant (either t_max or P_amb)

% Residual internal energy
finalGasEnergy.U_gas_residual = gasMass * Cv_gas * T_gas_final;

% Residual kinetic energy (estimate)
gasVelocityFactor_k = 1.7; % Same factor used before
v_gas_final_est = gasVelocityFactor_k * v_proj_final;
finalGasEnergy.KE_gas_estimated = 0.5 * gasMass * v_gas_final_est^2;

% Heat lost during expansion
finalGasEnergy.Q_loss_during_expansion = expansionResults.Q_loss_exp(end);

fprintf('--- Gas Expansion Results ---\n');
fprintf('  Final Gas Temp: %.1f K\n', T_gas_final);
fprintf('  Final Projectile Velocity: %.1f m/s\n', v_proj_final);
fprintf('  Final Gas Pressure: %.2f MPa\n', expansionResults.P_gas_calc(end)/1e6);
fprintf('  Residual Gas Internal Energy: %.1f kJ\n', finalGasEnergy.U_gas_residual / 1000);
fprintf('  Residual Gas Kinetic Energy (est. k=%.1f): %.1f kJ\n', gasVelocityFactor_k, finalGasEnergy.KE_gas_estimated / 1000);
fprintf('  Heat Lost During Expansion: %.1f kJ\n', finalGasEnergy.Q_loss_during_expansion / 1000);
fprintf('  Total Final Gas Energy (U+KE): %.1f kJ\n', (finalGasEnergy.U_gas_residual + finalGasEnergy.KE_gas_estimated)/1000);


end % End of function simulateGasExpansion


% =========================================================================
% --- MODIFIED Local ODE Function for Expansion Only ---
% =========================================================================
% MODIFIED Signature: Removed h_conv, added h_conv_base, h_conv_pref, h_conv_n
function dydt = odesExpansionOnly(t, y, m_gas, Cv_gas, m_p, A_b, D_bore, V0, R_gas, b_i, T_barrel, h_conv_base, h_conv_pref, h_conv_n)
    % State vector y = [T_gas; x_p; v_p; Q_loss_exp]

    T_gas = max(T_barrel, y(1)); % Gas temperature >= barrel temp
    x_p   = y(2);                % Allow any projectile position
    v_p   = max(0, y(3));        % Velocity >= 0

    % Intermediate Calculations
    V_inst = V0 + A_b * x_p;
    V_inst = max(1e-12, V_inst); % Instantaneous volume
    effective_volume = V_inst - m_gas * b_i;
    effective_volume = max(1e-9, effective_volume); % Effective volume (Noble-Able)

    P_gas = (m_gas * R_gas * T_gas) / effective_volume;
    P_gas = max(1e3, P_gas); % Minimum pressure (e.g., 1 kPa) to avoid issues

    % --- MODIFIED: Calculate Instantaneous h_conv ---
    % Use max(eps, ...) to prevent issues if P_gas is zero or negative briefly
    h_conv_current = h_conv_base * ( max(eps, P_gas / h_conv_pref) )^h_conv_n;
    h_conv_current = max(0, h_conv_current); % Ensure non-negative
    % --- End MODIFIED ---

    % Calculate Derivatives
    dydt = zeros(4, 1);

    % dT_gas/dt (only expansion work and heat loss)
    work_done_rate = P_gas * A_b * v_p; % P*dV/dt = P*A*v

    % Estimate heat transfer area (simple model)
    chamberAreaEstimate = 1e-3; % Use same estimate as main simulation (or make parameter)
    exposedCylinderArea = pi * D_bore * max(0, x_p); % Use max(0, x_p) for area
    contactArea = chamberAreaEstimate + exposedCylinderArea;

    % Calculate dQ_loss_dt using the dynamic h_conv_current
    dQ_loss_dt = 0; % Initialize
    if T_gas > T_barrel
        dQ_loss_dt = h_conv_current * contactArea * (T_gas - T_barrel); % USE h_conv_current
        dQ_loss_dt = max(0, dQ_loss_dt); % Heat loss rate cannot be negative
    end

    if v_p == 0
        work_done_rate = 0;
        % Keep dQ_loss_dt as calculated if T_gas > T_barrel
    end

    dT_dt_numerator = - work_done_rate - dQ_loss_dt; % dU = dQ - dW -> m*Cv*dT = -dQ_loss - P*dV
    dT_dt_denominator = m_gas * Cv_gas;

    if abs(dT_dt_denominator) < 1e-9
        dydt(1) = 0; % dT_gas/dt
    else
        dydt(1) = dT_dt_numerator / dT_dt_denominator;
    end

    % dx_p/dt
    dydt(2) = v_p; % dx_p/dt

    % dv_p/dt (only pressure force, friction neglected here)
    if m_p > 1e-9
        dydt(3) = (P_gas * A_b) / m_p; % dv_p/dt
    else
        dydt(3) = 0;
    end

    % dQ_loss_exp/dt
    dydt(4) = dQ_loss_dt; % dQ_loss/dt

end % End of function odesExpansionOnly


% =========================================================================
% --- Atmospheric Pressure Event Function (Unchanged) ---
% =========================================================================
function [value, isterminal, direction] = atmosphericPressureEvent(t, y, m_gas, R_gas, b_i, V0, A_b, P_amb)
% Event function to stop the simulation when P_gas <= P_ambient
    T_gas = y(1); x_p = y(2);
    V_inst = V0 + A_b * x_p; V_inst = max(1e-12, V_inst);
    effective_volume = V_inst - m_gas * b_i; effective_volume = max(1e-9, effective_volume);
    if T_gas < 0 || effective_volume <= 1e-10, P_gas_calc = P_amb; else P_gas_calc = (m_gas * R_gas * T_gas) / effective_volume; end
    P_gas_calc = max(0, P_gas_calc);
    value = P_gas_calc - P_amb;
    isterminal = 1; direction = -1;
end % End of function atmosphericPressureEvent