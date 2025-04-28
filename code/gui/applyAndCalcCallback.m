% gui/applyAndCalcCallback.m (or local function in ballisticSimulatorGUI.m)
% --- FINAL VERSION with External Physics Config Loading & Restored Error Handling ---
function success = applyAndCalcCallback(hFigure)
    success = false; % Default to failure
    handles = guidata(hFigure);
    % Use try-catch for logStatus in case log window handle is invalid during initial call
    try logStatus('"Apply & Calc" requested...'); catch, disp('"Apply & Calc" requested...'); end

    % Check prerequisites
    if ~isfield(handles,'isComponentsLoaded') || ~handles.isComponentsLoaded || ~isfield(handles,'componentData') || isempty(handles.componentData)
        try logStatus('Error: Load components before Apply & Calc.'); catch, disp('Error: Load components before Apply & Calc.'); end
        uiwait(msgbox('Load components first using the "Load Selected Components" button.', 'Load Required', 'warn'));
        return; % Exit if components not loaded
    end

    try % Main try block for the whole calculation process
        cd = handles.componentData; % Loaded component data (SI units)
        C = handles.CONV;           % Conversion constants struct
        inputs = struct();          % Struct to hold validated GUI inputs

        % --- Load Physics Configuration ---
        physicsConfigFile = fullfile('config', 'defaultPhysicsConfig.m'); % Adjust path if needed
        if ~exist(physicsConfigFile, 'file')
            error('Physics configuration file not found: %s', physicsConfigFile);
        end
        try
            clear physicsConfig; % Ensure clean workspace for the struct
            run(physicsConfigFile); % Execute the script to define 'physicsConfig'
            if ~exist('physicsConfig', 'var') || ~isstruct(physicsConfig)
                error('Script "%s" did not define the "physicsConfig" struct.', physicsConfigFile);
            end
             try logStatus('Loaded physics configuration from defaultPhysicsConfig.m'); catch, disp('Loaded physics configuration from defaultPhysicsConfig.m'); end
        catch ME_physicsLoad
            error('Error loading or running physics config file "%s": %s', physicsConfigFile, ME_physicsLoad.message);
        end

        % --- Define fields to read from GUI ---
        names = { 'BulletName', 'BulletMass_gr', 'BulletLength_in', 'BulletDiameter_in', ...
                  'CartridgeName', 'CaseLength_in', 'MaxCaseCapacity_grH2O', 'BoreDiameter_in', ...
                  'SeatingDepth_in', ...
                  'PropellantCharge_gr', 'BarrelLength_in', 'TwistRate_inPerTurn', ...
                  'InitialTemperature_K', 'AmbientPressure_Pa', 'ShotStartPressure_MPa' };

        valid = true; % Flag for overall input validity
        errorMsg = {}; % Cell array for validation error messages
        isNumericField = @(tag) ~any(strcmpi(tag, {'BulletName', 'CartridgeName'}));

        % --- Validate User Inputs from GUI ---
        try logStatus('Validating GUI inputs...'); catch, disp('Validating GUI inputs...'); end
        for i=1:length(names)
             tagSuffix = names{i};
             editHandleName = ['edit' tagSuffix];
             if ~isfield(handles, editHandleName) || ~ishandle(handles.(editHandleName)), error('GUI Handle not found for editable field: %s.', editHandleName); end
             h = handles.(editHandleName);
             valStr = get(h, 'String');
             if isNumericField(tagSuffix)
                 valNum = str2double(valStr); inputs.(tagSuffix) = valNum; isValidField = true; msg = '';
                 if isnan(valNum) || ~isreal(valNum), isValidField=false; msg='Invalid Number'; end
                 if isValidField % Specific validation checks
                     if any(strcmpi(tagSuffix,{'BulletMass_gr','BulletLength_in','BulletDiameter_in','CaseLength_in','MaxCaseCapacity_grH2O','BoreDiameter_in','SeatingDepth_in','BarrelLength_in'})) && valNum<=0, isValidField=false; msg='Must be > 0'; end
                     if strcmpi(tagSuffix,'PropellantCharge_gr') && valNum<0, isValidField=false; msg='Must be >= 0 gr'; end
                     if strcmpi(tagSuffix,'TwistRate_inPerTurn') && valNum==0, isValidField=false; msg='Cannot be 0 in/turn'; end
                     if strcmpi(tagSuffix,'InitialTemperature_K') && valNum<=0, isValidField=false; msg='Must be > 0 K'; end
                     if strcmpi(tagSuffix,'AmbientPressure_Pa') && valNum<0, isValidField=false; msg='Must be >= 0 Pa'; end
                     if strcmpi(tagSuffix,'ShotStartPressure_MPa') && valNum<=0, isValidField=false; msg='Must be > 0 MPa'; end
                 end
                 if ~isValidField, set(h, 'BackgroundColor', [1 0.8 0.8]); valid=false; errorMsg{end+1}=sprintf('%s: %s (%s)', tagSuffix, msg, valStr); else set(h, 'BackgroundColor', [1 1 1]); end
             else % String fields
                 cleanStr = strtrim(valStr); inputs.(tagSuffix) = cleanStr;
                 if isempty(cleanStr), set(h, 'BackgroundColor', [1 0.8 0.8]); valid=false; errorMsg{end+1}=sprintf('%s: Cannot be empty', tagSuffix); else set(h, 'BackgroundColor', [1 1 1]); end
             end
        end % End input validation loop

        if ~valid
            error('Invalid user inputs found: %s', strjoin(errorMsg,'; '));
        end
        if inputs.PropellantCharge_gr == 0, try logStatus('Warning: Propellant charge is zero.'); catch, disp('Warning: Propellant charge is zero.'); end; end

        % --- Perform Conversions to SI base units ---
        try logStatus('Converting inputs to SI units...'); catch, disp('Converting inputs to SI units...'); end
        propCharge_kg = inputs.PropellantCharge_gr * C.KG_PER_GRAIN;
        barrelLength_m = inputs.BarrelLength_in * C.METERS_PER_INCH;
        twistRate_inPerTurn = inputs.TwistRate_inPerTurn;
        initialTemp_K = inputs.InitialTemperature_K;
        ambientPress_Pa = inputs.AmbientPressure_Pa;
        shotStartPressure_Pa = inputs.ShotStartPressure_MPa * C.PA_PER_MPA;
        bulletMass_kg = inputs.BulletMass_gr * C.KG_PER_GRAIN;
        bulletDiameter_m = inputs.BulletDiameter_in * C.METERS_PER_INCH;
        bulletLength_m = inputs.BulletLength_in * C.METERS_PER_INCH;
        caseLength_m = inputs.CaseLength_in * C.METERS_PER_INCH;
        maxCapacity_m3 = inputs.MaxCaseCapacity_grH2O * C.M3_PER_GRH2O;
        boreDiameter_m = inputs.BoreDiameter_in * C.METERS_PER_INCH;
        seatingDepth_m = getfield_safe(cd.bullet, 'seatingDepth_m', NaN);
        if isnan(seatingDepth_m), error('Internal Error: Seating depth (seatingDepth_m) not found in loaded componentData.bullet struct.'); end

        % --- Calculate OAL (Overall Length) ---
        calculated_OAL_m = caseLength_m + bulletLength_m - seatingDepth_m;
        calculated_OAL_in = calculated_OAL_m * C.INCHES_PER_METER;
        try logStatus(sprintf('Using Seating Depth=%.3f in -> Calculated OAL=%.3f in.', inputs.SeatingDepth_in, calculated_OAL_in)); catch, disp(sprintf('Using Seating Depth=%.3f in -> Calculated OAL=%.3f in.', inputs.SeatingDepth_in, calculated_OAL_in)); end

        % --- Calculations Volume (SI) ---
        try logStatus('Calculating volumes...'); catch, disp('Calculating volumes...'); end
        powderDensity_kgm3 = getfield_safe(cd.powder, 'propellantDensity_rho_s', NaN);
        if isnan(powderDensity_kgm3) || powderDensity_kgm3 <= 0, error('Invalid powder density (propellantDensity_rho_s = %.3f) from loaded powder data.', powderDensity_kgm3); end
        powderVolume_m3 = 0; if propCharge_kg > 0, powderVolume_m3 = propCharge_kg / powderDensity_kgm3; end
        bulletArea_m2 = pi * (bulletDiameter_m / 2)^2;
        seatedBulletVolume_m3 = bulletArea_m2 * seatingDepth_m;
        qlUseableCapacity_m3 = maxCapacity_m3 - seatedBulletVolume_m3;
        initialFreeVolumeV0_m3 = qlUseableCapacity_m3 - powderVolume_m3;
        if qlUseableCapacity_m3 <= 0, warning('Calculated Useable Case Capacity (Max Capacity - Seated Bullet Volume) is %.3e m^3 (<= 0). Check inputs.', qlUseableCapacity_m3); end
        if initialFreeVolumeV0_m3 <= 1e-12
            if propCharge_kg > 0, error('Calculated Initial Free Volume V0 (Useable Capacity - Powder Volume) is %.3e m^3 (<= 0). Check inputs.', initialFreeVolumeV0_m3); else initialFreeVolumeV0_m3 = qlUseableCapacity_m3; try logStatus('V0 set to useable case capacity as propellant charge is zero.'); catch, disp('V0 set to useable case capacity as propellant charge is zero.'); end; end
        end
        bulletTravel_m = barrelLength_m - seatingDepth_m;
        if bulletTravel_m <= 0, warning('Calculated bullet travel (Barrel Length - Seating Depth) is %.3f m (<= 0). Check inputs.', bulletTravel_m); end

        % --- Update CALCULATED Read-Only Display Fields in GUI ---
        try logStatus('Updating calculated display fields...'); catch, disp('Updating calculated display fields...'); end
        set(handles.textVolumeOccupiedBullet_grH2O, 'String', sprintf('%.3f', seatedBulletVolume_m3 * C.GRH2O_PER_M3));
        set(handles.textQLUseableCaseCapacity_grH2O, 'String', sprintf('%.3f', qlUseableCapacity_m3 * C.GRH2O_PER_M3));
        set(handles.textInitialFreeVolumeV0_grH2O, 'String', sprintf('%.3f', initialFreeVolumeV0_m3 * C.GRH2O_PER_M3));
        set(handles.textCartridgeOAL_in, 'String', sprintf('%.3f', calculated_OAL_in));
        set(handles.textBulletTravel_in, 'String', sprintf('%.3f', bulletTravel_m * C.INCHES_PER_METER));

        % --- ASSEMBLE FINAL 'parameters' STRUCT for simulation ---
        try logStatus('Assembling final parameters struct for simulation...'); catch, disp('Assembling final parameters struct for simulation...'); end
        simParams = struct();
        simParams.loadedCartridgeName = inputs.CartridgeName;
        simParams.loadedBulletName = inputs.BulletName;
        powderList = get(handles.popupmenuPowder, 'String'); powderIndex = get(handles.popupmenuPowder, 'Value'); selectedPowderName = '';
        if iscell(powderList)
            if powderIndex > 0 && powderIndex <= length(powderList), selectedPowderName = powderList{powderIndex}; else error('Invalid powder selection index.'); end
        elseif ischar(powderList) && ~isempty(powderList) && ~contains(powderList, {'(','No ','Err'}), selectedPowderName = powderList; else error('No valid powder selected or available.'); end
        if isempty(selectedPowderName), error('Could not determine selected powder name.'); end
        simParams.loadedPowderName = selectedPowderName;
        simParams.projWeight_gr = inputs.BulletMass_gr;
        simParams.projDiameter_in = inputs.BulletDiameter_in;
        simParams.projMass_m = bulletMass_kg;
        simParams.projArea_Ab = bulletArea_m2;
        simParams.projMomentOfInertia_Ip = 0.5 * bulletMass_kg * (bulletDiameter_m / 2)^2; % Estimate
        simParams.initialPropellantMass_m = propCharge_kg;
        simParams.initialFreeVolume_V0 = initialFreeVolumeV0_m3;
        simParams.barrelLength_m = barrelLength_m;
        simParams.boreDiameter_Db = boreDiameter_m;
        simParams.twistRate_inPerTurn = twistRate_inPerTurn;
        if twistRate_inPerTurn ~= 0
            meters_per_turn = abs(twistRate_inPerTurn) * C.METERS_PER_INCH;
            simParams.twistRate_rad_m = sign(twistRate_inPerTurn) * (2 * pi) / meters_per_turn;
        else
            simParams.twistRate_rad_m = 0; try logStatus('Warn: Twist rate is zero, twistRate_rad_m set to 0.'); catch, disp('Warn: Twist rate is zero, twistRate_rad_m set to 0.'); end;
        end
        powderDataLoaded = cd.powder; powderFields = fieldnames(powderDataLoaded);
        try logStatus('Copying parameters from loaded powder data...'); catch, disp('Copying parameters from loaded powder data...'); end
        for k=1:length(powderFields), simParams.(powderFields{k}) = powderDataLoaded.(powderFields{k}); end
        R_univ = 8.314462;
        if isfield(simParams, 'molarMass_M_gas') && isnumeric(simParams.molarMass_M_gas) && simParams.molarMass_M_gas > 0
            simParams.specificGasConstant_R = R_univ / simParams.molarMass_M_gas;
             try logStatus(sprintf('Calculated specificGasConstant_R = %.2f J/(kg*K)', simParams.specificGasConstant_R)); catch, disp(sprintf('Calculated specificGasConstant_R = %.2f J/(kg*K)', simParams.specificGasConstant_R)); end
        else
            error('Molar mass (molarMass_M_gas=%.4e) missing or invalid in loaded powder data. Cannot calculate required specificGasConstant_R.', getfield_safe(simParams,'molarMass_M_gas',NaN));
        end
        sigma_powder = getfield_safe(powderDataLoaded,'specificSurfaceArea_sigma', NaN);
        if isnan(sigma_powder) || sigma_powder <= 0, error('Powder specific surface area (specificSurfaceArea_sigma=%.3f) invalid in loaded powder data.', sigma_powder); end
        simParams.specificSurfaceArea_sigma = sigma_powder;
        simParams.initialSurfaceArea_S0 = simParams.initialPropellantMass_m * simParams.specificSurfaceArea_sigma;
        if simParams.initialPropellantMass_m == 0, simParams.initialSurfaceArea_S0 = 0;
        elseif simParams.initialSurfaceArea_S0 <= 0, simParams.initialSurfaceArea_S0 = 1e-9; try logStatus('Warn: Calculated initialSurfaceArea_S0 was <= 0 despite non-zero charge, set to 1e-9 m^2.'); catch, disp('Warn: Calculated initialSurfaceArea_S0 was <= 0 despite non-zero charge, set to 1e-9 m^2.'); end; end
        simParams.bullet_bc_g1 = getfield_safe(cd.bullet, 'bc_g1', NaN);
        try logStatus(sprintf('Copied bullet_bc_g1 = %.3f from loaded bullet data.', simParams.bullet_bc_g1)); catch, disp(sprintf('Copied bullet_bc_g1 = %.3f from loaded bullet data.', simParams.bullet_bc_g1)); end
  % --- >>> ADD THIS LINE <<< ---
    simParams.formFactor_G1 = getfield_safe(cd.bullet, 'formFactor_g1', NaN); % Copy formFactor_G1
    try logStatus(sprintf('Copied formFactor_G1 = %.4f from loaded bullet data.', simParams.formFactor_G1)); catch, disp(sprintf('Copied formFactor_G1 = %.4f from loaded bullet data.', simParams.formFactor_G1)); end
    % --- >>> END ADDED LINE <<< ---
        % --- Assign Physics Params from Loaded Config ---
        try logStatus('Assigning physics parameters from loaded configuration...'); catch, disp('Assigning physics parameters from loaded configuration...'); end
        simParams.engravingTorque_Nm = getfield_safe(physicsConfig, 'engravingTorque_Nm', NaN);
        simParams.engravingEndPosition_m = getfield_safe(physicsConfig, 'engravingEndPosition_m', NaN);
        simParams.rotationalFrictionCoeff_muRot = getfield_safe(physicsConfig, 'rotationalFrictionCoeff_muRot', NaN);
        % REMOVE: simParams.convectiveHeatCoeff_h = getfield_safe(physicsConfig, 'convectiveHeatCoeff_h', NaN);
        % ADD THESE:
        simParams.h_conv_base_coefficient = getfield_safe(physicsConfig, 'h_conv_base_coefficient', NaN);
        simParams.h_conv_pressure_reference_Pa = getfield_safe(physicsConfig, 'h_conv_pressure_reference_Pa', NaN);
        simParams.h_conv_pressure_exponent = getfield_safe(physicsConfig, 'h_conv_pressure_exponent', NaN);
        % END ADD

        simParams.maxSafetyTime_s = getfield_safe(physicsConfig, 'maxSafetyTime_s', 0.004);
        simParams.boreResistance_travel_m = getfield_safe(physicsConfig, 'boreResistance_travel_m', []);
        simParams.boreResistance_pressure_Pa = getfield_safe(physicsConfig, 'boreResistance_pressure_Pa', []);

        % Validate loaded physics parameters (Update this check)
        if isnan(simParams.engravingTorque_Nm) || isnan(simParams.engravingEndPosition_m) || ...
           isnan(simParams.rotationalFrictionCoeff_muRot) || ...
           isnan(simParams.h_conv_base_coefficient) || isnan(simParams.h_conv_pressure_reference_Pa) || isnan(simParams.h_conv_pressure_exponent) || ... % Added h_conv params check
           simParams.h_conv_pressure_reference_Pa <= 0 || ... % Ensure reference pressure is positive
           isempty(simParams.boreResistance_travel_m) || isempty(simParams.boreResistance_pressure_Pa) || ...
           length(simParams.boreResistance_travel_m) ~= length(simParams.boreResistance_pressure_Pa)
             warning('One or more physics/bore resistance parameters loaded from config are invalid or missing. Check defaultPhysicsConfig.m.');
             error('Failed to load valid physics parameters from config file.');
        end

        % --- Add initial conditions and environmental params from inputs ---
        simParams.initialTemperature_K = initialTemp_K;
        simParams.ambientPressure = ambientPress_Pa;
        simParams.shotStartPressure_Pa = shotStartPressure_Pa;

        % --- Final Validation ---
  % required_for_odes = { ...
  %           'initialPropellantMass_m', 'initialTemperature_K', 'projMass_m', 'projArea_Ab', 'initialFreeVolume_V0', ...
  %           'propellantDensity_rho_s', 'covolume_b', 'impetus_F', ...
  %           'specificHeatRatio_gamma', 'specificGasConstant_R', ...
  %           'burnRateCoeff_a', 'burnRateExponent_beta', 'initialSurfaceArea_S0', 'formFunctionParam_theta', ...
  %           'projMomentOfInertia_Ip', 'twistRate_rad_m', ...
  %           'shotStartPressure_Pa', ...
  %           'boreResistance_travel_m', 'boreResistance_pressure_Pa', ...
  %           'engravingTorque_Nm', 'engravingEndPosition_m', ...
  %           'rotationalFrictionCoeff_muRot', ...
  %           % 'convectiveHeatCoeff_h', ... % <<< REMOVE THIS LINE <<<
  %           'h_conv_base_coefficient', 'h_conv_pressure_reference_Pa', 'h_conv_pressure_exponent', ... % Ensure these are present
  %           'boreDiameter_Db', ...
  %           'barrelLength_m', 'maxSafetyTime_s', 'ambientPressure' ...
  %           }; % Make sure this list exactly matches the parameters needed by simulationOdes.m
  % 
  %       missing_odes_fields = required_for_odes(~isfield(simParams, required_for_odes));
  %       if ~isempty(missing_odes_fields)
  %            error('Internal Error: Assembled simParams struct is missing fields required by simulation_odes: %s', strjoin(missing_odes_fields,', '));
  %       end
        try logStatus('Final parameter check passed. Struct ready for simulation.'); catch, disp('Final parameter check passed. Struct ready for simulation.'); end

        % --- Store final parameters struct in main GUI handles ---
        handles.parametri = simParams;
        handles.isSimRun = false;
        guidata(hFigure, handles); % Save updated handles

        try logStatus('Parameters applied & calculated successfully. Ready to run simulation.'); catch, disp('Parameters applied & calculated successfully. Ready to run simulation.'); end
        success = true; % Indicate success

    catch ME % --- RESTORED Original CATCH block ---
        % Clear params on error to prevent running with bad data
        handles.parametri = [];
        handles.isSimRun = false;
        try guidata(hFigure, handles); catch, end % Try to save cleared state

        % Log and display error using original, more detailed methods
        msgError = ['ERROR during Apply & Calc: ', ME.message];
        try logStatus(msgError); catch, disp(msgError); end % Log to window or console
        fprintf(2,'Stack trace for Apply & Calc error in %s (Line %d):\n', ME.stack(1).file, ME.stack(1).line); % Print stack trace to command window (stderr)
        disp(ME.getReport('basic')); % Display basic error report in command window (stdout)
        uiwait(msgbox(sprintf('Error during Apply & Calculation:\n%s\n(Check log/console for details)', ME.message), 'Calculation Error', 'error')); % Show message box

        % Attempt to reset calculated display fields in GUI
        try
             set(handles.textVolumeOccupiedBullet_grH2O, 'String', 'ERROR');
             set(handles.textQLUseableCaseCapacity_grH2O, 'String', 'ERROR');
             set(handles.textInitialFreeVolumeV0_grH2O, 'String', 'ERROR');
             set(handles.textCartridgeOAL_in, 'String', 'ERROR');
             set(handles.textBulletTravel_in, 'String', 'ERROR');
        catch ME_gui_reset
            msgResetError = ['Error resetting calculated fields after Apply&Calc error: ', ME_gui_reset.message];
            try logStatus(msgResetError); catch, disp(msgResetError); end
        end
        success = false; % Indicate failure
    end % End main try-catch

end % End function applyAndCalcCallback

% --- Local helper function ---
function val = getfield_safe(s, field, default)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field)), val = s.(field); else val = default; end
end