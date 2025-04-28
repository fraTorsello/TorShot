    % =========================================================
    % Internal Ballistics Simulator - Refactored Version
    % =========================================================
    % core/simulationOdes.m
    % Function defining the system of Ordinary Differential Equations (ODEs)
    % for the internal ballistics model.
    % Calculates the time derivatives of the state variables.
    % Includes: propellant combustion, projectile motion (linear & rotational),
    %           gas thermodynamics (Noble-Able), tabulated bore resistance,
    %           simple heat loss, improved energy balance.
    %
    % NOTES:
    % - USES TABULATED BORE RESISTANCE br(x_p) instead of simple mu_friction.
    % - Implements Lagrange Gradient for base pressure P_base.
    % - dT/dt balance includes d(KE_gas)/dt and d(KE_rot)/dt terms.
    % =========================================================
    
    function dydt = simulationOdes(t, y, params)
    % INPUTS:
    %   t: Current time (scalar)
    %   y: State vector [m_prop_rem; T_gas; x_p; v_p; omega; W_bore_res; Q_loss]
    %      y(1): m_prop_rem         - Remaining propellant mass [kg]
    %      y(2): T_gas              - Average gas temperature [K]
    %      y(3): x_p                - Projectile position [m]
    %      y(4): v_p                - Projectile linear velocity [m/s]
    %      y(5): omega              - Projectile angular velocity [rad/s]
    %      y(6): W_bore_res         - Cumulative work against bore resistance [J]
    %      y(7): Q_loss             - Cumulative heat lost to walls [J]
    %   params: Struct containing all simulation parameters.
    %           MUST include params.boreResistance_travel_m and params.boreResistance_pressure_Pa
    % OUTPUT:
    %   dydt: Vector of derivatives [dm_prop_dt; dT_dt; dx_dt; dv_dt; domega_dt; dW_br_dt; dQ_loss_dt]
    
    % --- Extract State Variables from Vector y ---
    propellantMassRemaining = y(1);
    gasTemperature = y(2);
    projPosition = y(3);
    projVelocity = y(4);
    projAngularVelocity = y(5);
    workBoreResistance = y(6);
    heatLossCumulative = y(7);
    
    % --- Extract Necessary Parameters from params Struct ---
    % Main parameters
    initialPropellantMass = params.initialPropellantMass_m; % Used as charge mass C_T
    projMass_m = params.projMass_m;
    projArea_Ab = params.projArea_Ab;
    initialFreeVolume_V0 = params.initialFreeVolume_V0;
    propellantDensity_rho_s = params.propellantDensity_rho_s;
    covolume_b = params.covolume_b;
    specificHeatRatio_gamma = params.specificHeatRatio_gamma;
    specificGasConstant_R = params.specificGasConstant_R;
    impetus_F = params.impetus_F;
    % Combustion
    burnRateCoeff_a = params.burnRateCoeff_a;
    burnRateExponent_beta = params.burnRateExponent_beta;
    initialSurfaceArea_S0 = params.initialSurfaceArea_S0;
    formFunctionParam_theta = params.formFunctionParam_theta;
    % Rotational motion
    projMomentOfInertia_Ip = params.projMomentOfInertia_Ip;
    twistRate_rad_m = params.twistRate_rad_m;
    engravingTorque_Nm = params.engravingTorque_Nm; % Could be part of bore resistance?
    engravingEndPosition_m = params.engravingEndPosition_m; % Could be related to bore resistance?
    rotationalFrictionCoeff_muRot = params.rotationalFrictionCoeff_muRot; % Kept separate for now
    boreDiameter_Db = params.boreDiameter_Db;
    % Initial conditions and Resistance
    shotStartPressure_Pa = params.shotStartPressure_Pa;
    boreResistance_travel_m = params.boreResistance_travel_m; % Vector of positions for br [m]
    boreResistance_pressure_Pa = params.boreResistance_pressure_Pa; % Vector of resisting pressures br [Pa]
    % Ambient pressure handling
    if isfield(params, 'ambientPressure'), ambientPressure = params.ambientPressure; else, warning('ID:ambientP','ambientPressure not found in params, using 0 Pa.'); ambientPressure = 0.0; end
    if ~isnumeric(ambientPressure) || ~isscalar(ambientPressure) || ambientPressure < 0, warning('ID:ambientPInv','ambientPressure invalid, using 0 Pa.'); ambientPressure = 0.0; end
    % Heat Loss parameters
    h_conv_base = params.h_conv_base_coefficient;
    h_conv_pref = params.h_conv_pressure_reference_Pa;
    h_conv_n    = params.h_conv_pressure_exponent;
    initialBarrelTemp_K = params.initialTemperature_K;
    
    % --- Basic Input Validation for Bore Resistance (br) ---
    if length(boreResistance_travel_m) ~= length(boreResistance_pressure_Pa) || length(boreResistance_travel_m) < 2
        error('Vectors boreResistance_travel_m and boreResistance_pressure_Pa must have the same length (>= 2).');
    end
    if ~issorted(boreResistance_travel_m) || any(boreResistance_travel_m < 0)
         error('Vector boreResistance_travel_m must be monotonically increasing and non-negative.');
    end
    % Optional: check if bore resistance pressure is always non-negative
    % if any(boreResistance_pressure_Pa < 0)
    %     warning('ID:brNeg','Vector boreResistance_pressure_Pa contains negative values. They will be limited to 0.');
    %     boreResistance_pressure_Pa(boreResistance_pressure_Pa < 0) = 0;
    % end
    
    % --- Protections for Non-Physical Values ---
    propellantMassRemaining = max(1e-12, propellantMassRemaining);
    gasTemperature = max(initialBarrelTemp_K, gasTemperature);
    projVelocity = max(0, projVelocity);
    projAngularVelocity = max(0, projAngularVelocity);
    % Ensure projPosition is within or at the start of the defined range for br
    projPosition = max(boreResistance_travel_m(1), projPosition);
    
    % --- Fundamental Intermediate Calculations ---
    gasMass = initialPropellantMass - propellantMassRemaining; % Current gas mass [kg]
    gasMass = max(1e-12, gasMass); % Avoid zero mass
    
    instantaneousVolume = initialFreeVolume_V0 + projArea_Ab * projPosition; % Available volume [m^3]
    instantaneousVolume = max(1e-12, instantaneousVolume);
    
    % Specific heat at constant volume (Cv)
    if specificHeatRatio_gamma <= 1, error('specificHeatRatio_gamma must be > 1 (value: %.3f)', specificHeatRatio_gamma); end
    Cv_gas = specificGasConstant_R / (specificHeatRatio_gamma - 1); % [J/(kg*K)]
    Cv_gas = max(1e-6, Cv_gas); % Avoid zero or negative Cv
    
    % Average Gas Pressure (Noble-Able Law)
    effectiveVolume = instantaneousVolume - gasMass * covolume_b;
    effectiveVolume = max(1e-9, effectiveVolume); % Avoid non-physical volume
    avgGasPressure = (gasMass * specificGasConstant_R * gasTemperature) / effectiveVolume; % [Pa]
    avgGasPressure = max(1e3, avgGasPressure); % Avoid negative or extremely low pressures
    
    % --- Calculate Propellant Combustion Rate ---
    burningSurfaceArea = 0;
    burnRate = 0;
    dm_prop_dt = 0; % Propellant mass derivative [kg/s]
    dm_gas_dt = 0;  % Gas mass derivative [kg/s]
    
    if propellantMassRemaining > 1e-9 * initialPropellantMass % Only burn if significant propellant remains
        massFractionBurned_f = (initialPropellantMass - propellantMassRemaining) / initialPropellantMass;
        massFractionBurned_f = max(0, min(massFractionBurned_f, 1)); % Clamp between 0 and 1
        formFunction_Zf = (1 - massFractionBurned_f) * (1 + formFunctionParam_theta * massFractionBurned_f);
        formFunction_Zf = max(0, formFunction_Zf); % Ensure non-negative surface area factor
        burningSurfaceArea = initialSurfaceArea_S0 * formFunction_Zf; % [m^2]
    
        burnRate = burnRateCoeff_a * (avgGasPressure^burnRateExponent_beta); % Linear burn rate [m/s]
        burnRate = max(0, burnRate);
    
        dm_prop_dt = -propellantDensity_rho_s * burningSurfaceArea * burnRate; % [kg/s]
        dm_prop_dt = min(0, dm_prop_dt); % Derivative is negative or zero
        dm_gas_dt = -dm_prop_dt; % Mass conservation
    end
    
    % --- Calculate Resistances, Pressures, Torques ---
    
    % Interpolate Tabulated Bore Resistance (br)
    % 'extrap' uses the value at the nearest endpoint if projPosition is outside the range
    % Ensure vectors are columns for interp1 if needed
    currentBoreResistance_Pa = interp1(boreResistance_travel_m(:), boreResistance_pressure_Pa(:), projPosition, 'linear', 'extrap');
    currentBoreResistance_Pa = max(0, currentBoreResistance_Pa); % Force br >= 0 Pa
    
    % External Pressure 'Pg' (ambient pressure in front of projectile)
    externalPressure_Pg = ambientPressure; % [Pa]
    
    % Projectile Base Pressure (P_base) - Lagrange Gradient
    chargeMass_CT = initialPropellantMass; % Approximate total charge mass C_T with initial propellant mass
    lagrangeNumerator = avgGasPressure + (chargeMass_CT * (currentBoreResistance_Pa + externalPressure_Pg)) / (3.0 * projMass_m);
    lagrangeDenominator = 1.0 + (chargeMass_CT / (3.0 * projMass_m));
    if abs(lagrangeDenominator) < 1e-9
        projBasePressure = avgGasPressure; % Fallback if denominator is near zero
    else
        projBasePressure = lagrangeNumerator / lagrangeDenominator;
    end
    projBasePressure = max(0, projBasePressure); % [Pa]
    
    % Torques acting on the projectile
    % Note: Still using avgGasPressure for torques for simplicity/previous consistency.
    % Rotational friction is still tied to mu_rot. If br includes everything, remove mu_rot/torqueRotFriction.
    avgPressureForce = avgGasPressure * projArea_Ab;
    boreRadius_r = boreDiameter_Db / 2.0;
    riflingTorque = avgPressureForce * twistRate_rad_m * boreRadius_r^2;
    engravingTorque = (projPosition < engravingEndPosition_m) * engravingTorque_Nm; % May be redundant if included in br(x)
    rotationalFrictionTorque = avgPressureForce * rotationalFrictionCoeff_muRot * boreRadius_r; % Separate rotational friction
    resistiveTorque = max(0, engravingTorque + rotationalFrictionTorque); % Consider removing engravingTorque if part of br
    netTorque = riflingTorque - resistiveTorque;
    
    % --- Calculate Derivatives of State Variables ---
    dydt = zeros(7, 1);
    
    % 1. Propellant Mass Derivative [kg/s]
    dydt(1) = dm_prop_dt;
    
    % 2. Projectile Position Derivative (Velocity) [m/s]
    dx_dt = projVelocity;
    dydt(3) = dx_dt;
    
    % 3. Projectile Linear Velocity Derivative (Acceleration) [m/s^2]
    dv_dt = 0;
    % Net Accelerating Force [N] = Base Pressure Force - Bore Resistance Force - External Pressure Force
    netLinearForce = (projBasePressure * projArea_Ab) - (currentBoreResistance_Pa * projArea_Ab) - (externalPressure_Pg * projArea_Ab);
    % Conditions for acceleration (avgGasPressure > shotStartPressure only for initial start?)
    if netLinearForce > 0 && avgGasPressure > shotStartPressure_Pa % Starts moving only above shotStartPressure
       if projMass_m > 1e-9, dv_dt = netLinearForce / projMass_m; end
    elseif netLinearForce < 0 && projVelocity > 0 % Allows deceleration
        if projMass_m > 1e-9, dv_dt = netLinearForce / projMass_m; end
    else % Stationary or force insufficient/negative from rest
        dv_dt = 0;
    end
    dydt(4) = dv_dt;
    
    % 4. Projectile Angular Velocity Derivative (Angular Acceleration) [rad/s^2]
    domega_dt = 0;
    if netTorque > 0 && avgGasPressure > shotStartPressure_Pa % Accelerates only if net torque positive and pressure sufficient
        if projMomentOfInertia_Ip > 1e-12, domega_dt = netTorque / projMomentOfInertia_Ip; end
    elseif netTorque < 0 && projAngularVelocity > 0 % Allows angular deceleration
        if projMomentOfInertia_Ip > 1e-12, domega_dt = netTorque / projMomentOfInertia_Ip; end
    else
        domega_dt = 0;
    end
    dydt(5) = domega_dt;
    
    % 5. Derivative of Work against Bore Resistance (dW_br/dt) [J/s = W]
    if projVelocity > 0
        boreResistanceForce = currentBoreResistance_Pa * projArea_Ab; % Instantaneous resisting force [N]
        dW_br_dt = boreResistanceForce * projVelocity; % Power dissipated [W]
    else
        dW_br_dt = 0; % No work if stationary
    end
    dydt(6) = dW_br_dt;
    
    % 6. Gas Temperature Derivative (REVISED Energy Balance) [K/s]
    dT_dt = 0;
    dQ_loss_dt = 0; % Heat loss rate to walls
    if gasMass > 1e-9 && Cv_gas > 0
        % Energy rate terms (Watts)
        % Energy input from combustion: d(m_gas)/dt * (Internal Energy per mass - Cv*T)
        energyInputRate = dm_gas_dt * (impetus_F - Cv_gas * gasTemperature);
        % Work done rate by gas on projectile: P_base * A_b * v_p
        workDoneRate = projBasePressure * projArea_Ab * projVelocity;
        workDoneRate = max(0, workDoneRate); % Work done by gas is positive or zero
    
        % Heat loss rate (Simple model)
        chamberAreaEstimate = 1e-3; % Estimate fixed chamber area - improve if needed
        exposedCylinderArea = pi * boreDiameter_Db * max(0, projPosition); % Area of barrel exposed by projectile
        contactArea = chamberAreaEstimate + exposedCylinderArea;
        % Calculate Instantaneous h_conv
h_conv_current = h_conv_base * ( max(eps, avgGasPressure / h_conv_pref) )^h_conv_n;
h_conv_current = max(0, h_conv_current); % Ensure non-negative
%disp(['🔥 Heat transfer: H used = ', num2str(h_conv_current), ' W/m²·K']);
        if gasTemperature > initialBarrelTemp_K
          dQ_loss_dt = h_conv_current * contactArea * (gasTemperature - initialBarrelTemp_K);
        dQ_loss_dt = max(0, dQ_loss_dt); % Heat loss rate >= 0
        else
            dQ_loss_dt = 0;
        end
    
        % Rate of change of Gas Kinetic Energy (Lagrange Approximation)
        % d(KE_gas)/dt = d(1/2 * (C_T/3) * v_p^2)/dt = (C_T/3) * v_p * dv_dt
        dKE_gas_dt_lagrange = (chargeMass_CT / 3.0) * projVelocity * dv_dt; % Requires dv_dt calculated above
    
        % Rate of change of Projectile Rotational Kinetic Energy
        % d(KE_rot)/dt = d(1/2 * I_p * omega^2)/dt = I_p * omega * domega_dt
        dKE_rot_dt = projMomentOfInertia_Ip * projAngularVelocity * domega_dt; % Requires domega_dt calculated above
    
        % Energy Balance: dU/dt = dQ_in/dt - dW_out/dt - dQ_loss/dt - dKE_gas/dt - dKE_rot/dt
        % m*Cv*dT/dt = energyInputRate - workDoneRate - dQ_loss_dt - dKE_gas_dt_lagrange - dKE_rot_dt
        dT_dt_numerator = energyInputRate - workDoneRate - dQ_loss_dt - dKE_gas_dt_lagrange - dKE_rot_dt;
        dT_dt_denominator = gasMass * Cv_gas;
    
        if abs(dT_dt_denominator) < 1e-9
            dT_dt = 0;
        else
            dT_dt = dT_dt_numerator / dT_dt_denominator;
        end
    else
        dT_dt = 0;
        dQ_loss_dt = 0;
    end
    dydt(2) = dT_dt;
    
    % 7. Derivative of Cumulative Heat Lost (to walls) [J/s = W]
    dydt(7) = dQ_loss_dt;
    
    % --- End Derivative Calculation ---
    end % End function simulationOdes