% config/defaultPhysicsConfig.m
% Configuration file for physics model estimates and bore resistance

disp('Executing defaultPhysicsConfig.m');

physicsConfig = struct();

% --- Physics Model Estimates ---
% NOTE: These are general estimates and might need tuning per caliber/bullet type
physicsConfig.engravingTorque_Nm           = 20;     % Engraving torque [N*m] (Estimate)
physicsConfig.engravingEndPosition_m       = 0.05;    % Position where engraving is complete [m] (Estimate)
physicsConfig.rotationalFrictionCoeff_muRot= 0.0279;   % Rotational friction coefficient [-] (Estimate)
% physicsConfig.convectiveHeatCoeff_h        =0;   % Convective heat transfer coefficient [W/(m^2*K)] (Estimate)

physicsConfig.h_conv_base_coefficient    = 15e4;   % Base coefficient [W/(m^2*K)] (Example value)
physicsConfig.h_conv_pressure_reference_Pa = 100e6; % Reference pressure [Pa] (Example: 100 MPa)
physicsConfig.h_conv_pressure_exponent   = 0.8;   

% --- Bore Resistance Tab   le ---
% Defines the pressure opposing the bullet's motion along the barrel.
% Should ideally be specific to the bullet/barrel combination.
% Format: [travel_vector_m], [pressure_vector_Pa]
physicsConfig.boreResistance_travel_m    = [0.0, 0.002, 0.010, 0.050, 0.600]; % Travel distance [m]
physicsConfig.boreResistance_pressure_Pa = [0, 20e6, 10e6, 5e6,   5e6];   % Corresponding resistance pressure [Pa] (Example values!)

% --- Simulation Control (Optional - could remain in applyAndCalcCallback) ---
physicsConfig.maxSafetyTime_s            = 0.0040;  % Max simulation time [s]
% Note: shotStartPressure_MPa is currently a GUI input, so it's not included here.

disp('Struct "physicsConfig" defined.');

% DO NOT add 'clearvars' or 'clear' at the end