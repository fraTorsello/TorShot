% File: analyzeAndSaveResults.m
% Description: Postprocessing of simulation results. Saves results and parameters.
% MODIFIED: Calls calculateBarrelStresses and saves stress results.
% =========================================================

% MODIFIED: Removed 'componentData' input argument.
% MODIFIED: Added logic to load componentData just before TXT writing and MAT saving.
% MODIFIED: Calls calculateBarrelStresses and saves the output.
function analysisOutput = analyzeAndSaveResults(simulationResults, parameters) % Removed barrelParameters input
    disp('Starting results analysis and saving...');

    % === Input validation ===
    if ~isstruct(simulationResults) || ~isfield(simulationResults, 'timeS') || isempty(simulationResults.timeS)
        error('analyzeAndSaveResults: Invalid or empty simulationResults struct provided.');
    end
    % Check parameters struct minimally and for required names for reloading
    required_names = {'loadedCartridgeName', 'loadedBulletName', 'loadedPowderName'};
    if ~isstruct(parameters) || ~all(isfield(parameters, required_names))
         error('analyzeAndSaveResults: Invalid or incomplete parameters struct provided (missing component names).');
    end
    % Extract names needed later for potential reloading
    cartridgeName = parameters.loadedCartridgeName;
    bulletName = parameters.loadedBulletName;
    powderName = parameters.loadedPowderName;
    if isempty(cartridgeName) || isempty(bulletName) || isempty(powderName)
        warning('analyzeAndSaveResults: Component names in parameters struct are empty. Cannot load component data.');
        canLoadComponentData = false;
    else
        canLoadComponentData = true; % Assume we can try loading
    end

    % === Load Barrel Data (needed for stress calc and saving) ===
    barrelData = [];
    barrelDataLoaded = false;
    if canLoadComponentData && exist('loadBarrelData', 'file') == 2
        try
            selectedBarrelName = parameters.selectedBarrelName; % Assumes this was stored in parameters
            if isempty(selectedBarrelName)
                 warning('analyzeAndSaveResults: No selectedBarrelName found in parameters. Cannot load barrel data.');
            else
                 barrelData = loadBarrelData(selectedBarrelName);
                 if isstruct(barrelData) && isfield(barrelData, 'yieldStrength_Pa') % Basic check
                     barrelDataLoaded = true;
                     disp('Barrel data loaded for stress calculation and saving.');
                 else
                     warning('analyzeAndSaveResults: Loaded barrelData appears invalid.');
                 end
            end
        catch ME_load_barrel
             warning('analyzeAndSaveResults: Failed to load barrelData: %s.', ME_load_barrel.message);
        end
    end
    if ~barrelDataLoaded
         warning('analyzeAndSaveResults: Barrel data could not be loaded. Skipping stress calculation and saving barrel info.');
    end

    % === Extract data ===
    % (Keep existing extraction code)
    timeS = simulationResults.timeS;
    gasPressurePa = simulationResults.gasPressurePa;
    projectilePositionM = simulationResults.projectilePositionM;
    projectileVelocityMps = simulationResults.projectileVelocityMps;
    angularVelocityRadps = simulationResults.angularVelocityRadps;
    remainingPropellantMassKg = simulationResults.remainingPropellantMassKg;
    frictionWorkJ = simulationResults.frictionWorkJ;
    heatLossJ = simulationResults.heatLossJ;
    gasMassKg = simulationResults.gasMassKg;
    gasTemperatureK = simulationResults.gasTemperatureK;
    riflingTorqueNm = simulationResults.riflingTorqueNm;
    engravingTorqueNm = simulationResults.engravingTorqueNm;
    frictionTorqueNm = simulationResults.frictionTorqueNm;
    netTorqueNm = simulationResults.netTorqueNm;

    % === Calculate Stresses ===
    stressResults = []; % Initialize
    if barrelDataLoaded && exist('calculateBarrelStresses', 'file') == 2
        try
            stressResults = calculateBarrelStresses(simulationResults, barrelData);
        catch ME_stress
            warning('analyzeAndSaveResults: Error calculating barrel stresses: %s', ME_stress.message);
            stressResults = []; % Ensure it's empty on error
        end
    elseif barrelDataLoaded
         warning('analyzeAndSaveResults: calculateBarrelStresses.m function not found. Skipping stress calculation.');
    end

    % === Final state ===
    % (Keep existing code)
    endIndex = length(timeS);
    finalVelocitySim = projectileVelocityMps(endIndex);
    finalOmegaSimRadps = angularVelocityRadps(endIndex);
    finalOmegaSimRpm = finalOmegaSimRadps * 60 / (2*pi);
    maxPressureSimPa = max(gasPressurePa);
    exitTimeSimS = timeS(endIndex);
    finalPositionM = projectilePositionM(endIndex);

    % === Parameters ===
    % (Keep existing code)
    initialPropMass = parameters.initialPropellantMass_m;
    barrelLen = parameters.barrelLength_m;
    projMass = parameters.projMass_m;
    projInertia = parameters.projMomentOfInertia_Ip;
    impetus = parameters.impetus_F;
    gasConst = parameters.specificGasConstant_R;
    gamma = parameters.specificHeatRatio_gamma;
    caseNameToDisplay = parameters.loadedCartridgeName;

    propellantConsumedKg = max(0, initialPropMass - remainingPropellantMassKg(endIndex));
    totalHeatLossJ = heatLossJ(endIndex);
    totalFrictionWorkJ = frictionWorkJ(endIndex);

    % === Prepare text summary ===
    summaryLines = {};
    function addLine(str), summaryLines{end+1, 1} = str; end

    % --- Populate Simulation Summary Lines ---
    % (Keep existing code)
    addLine(sprintf('--- TORSHOT Simulation Summary (%s) ---', caseNameToDisplay));
    addLine(sprintf('Simulation Time:               %.4f ms', exitTimeSimS * 1000));
    addLine(sprintf('Final Position:                %.4f m (Barrel: %.4f m)', finalPositionM, barrelLen));
    if abs(finalPositionM - barrelLen) < 1e-3, addLine('  -> Barrel exit occurred.'); else addLine('  -> WARNING: Barrel exit NOT reached.'); end
    addLine(sprintf('Final Linear Velocity:         %.2f m/s', finalVelocitySim));
    addLine(sprintf('Final Angular Velocity:        %.0f RPM (%.2f rad/s)', finalOmegaSimRpm, finalOmegaSimRadps));
    addLine(sprintf('Maximum Pressure:              %.2f MPa', maxPressureSimPa / 1e6));
    if initialPropMass > 1e-12, addLine(sprintf('Propellant Consumed:           %.2f %%', (propellantConsumedKg / initialPropMass) * 100)); else addLine('Propellant Consumed:           N/A (Initial Mass = 0)'); end
    addLine(sprintf('Total Bore Resistance Work:    %.2f kJ', totalFrictionWorkJ / 1000));
    addLine(sprintf('Total Heat Loss to Barrel:     %.2f kJ', totalHeatLossJ / 1000));

    % --- Add Stress Summary Lines ---
    if ~isempty(stressResults)
         addLine(sprintf('Maximum Equivalent Stress:     %.2f MPa', stressResults.max_eq_stress_Pa / 1e6));
         addLine(sprintf('Minimum Safety Factor:         %.2f (at %.4f ms)', stressResults.min_safety_factor, stressResults.time_at_min_sf_s * 1000));
         addLine(sprintf('Barrel Yield Strength:         %.1f MPa', stressResults.yield_strength_Pa / 1e6));
    else
         addLine('Maximum Equivalent Stress:     N/A');
         addLine('Minimum Safety Factor:         N/A');
    end

    % === Energy balance Calculation ===
    % (Keep existing code)
    releasedCombustionEnergyJ = propellantConsumedKg * impetus;
    finalLinearKeJ = 0.5 * projMass * finalVelocitySim^2;
    finalRotationalKeJ = 0.5 * projInertia * finalOmegaSimRadps^2;
    finalGasInternalEnergyJ = 0;
    if ~isnan(gamma) && (gamma - 1) > 1e-9 && ~isnan(gasConst) && ~isnan(gasMassKg(endIndex)) && ~isnan(gasTemperatureK(endIndex))
        Cv = gasConst / (gamma - 1);
        finalGasInternalEnergyJ = gasMassKg(endIndex) * Cv * gasTemperatureK(endIndex);
    else warning('analyzeAndSaveResults: Could not calculate finalGasInternalEnergyJ.'); end
    energyValues = [releasedCombustionEnergyJ, totalFrictionWorkJ, finalLinearKeJ, finalRotationalKeJ, totalHeatLossJ, finalGasInternalEnergyJ];
    if any(isnan(energyValues)), unaccountedEnergyJ = NaN; warning('analyzeAndSaveResults: Unaccounted energy cannot be calculated accurately.'); else unaccountedEnergyJ = releasedCombustionEnergyJ - totalFrictionWorkJ - finalLinearKeJ - finalRotationalKeJ - totalHeatLossJ - finalGasInternalEnergyJ; end

    % --- Populate Energy Balance Summary Lines ---
    % (Keep existing code)
    addLine('');
    addLine('--- Approximate Energy Balance ---');
    if ~isnan(releasedCombustionEnergyJ) && releasedCombustionEnergyJ > 1e-9
        addLine(sprintf('Released Chemical Energy:      %.2f kJ', releasedCombustionEnergyJ / 1000));
        addLine(sprintf('  -> Linear Kinetic Energy:     %.2f kJ (%.1f%%)', finalLinearKeJ / 1000, (finalLinearKeJ / releasedCombustionEnergyJ) * 100));
        addLine(sprintf('  -> Rotational Kinetic Energy: %.3f kJ (%.2f%%)', finalRotationalKeJ / 1000, (finalRotationalKeJ / releasedCombustionEnergyJ) * 100));
        addLine(sprintf('  -> Bore Resistance Work Loss: %.2f kJ (%.1f%%)', totalFrictionWorkJ / 1000, (totalFrictionWorkJ / releasedCombustionEnergyJ) * 100));
        addLine(sprintf('  -> Heat Loss to Barrel:       %.2f kJ (%.1f%%)', totalHeatLossJ / 1000, (totalHeatLossJ / releasedCombustionEnergyJ) * 100));
        addLine(sprintf('  -> Residual Gas Energy:       %.2f kJ (%.1f%%)', finalGasInternalEnergyJ / 1000, (finalGasInternalEnergyJ / releasedCombustionEnergyJ) * 100));
        addLine(sprintf('  -> Unaccounted Energy:        %.3f kJ (%.2f%%)', unaccountedEnergyJ / 1000, (unaccountedEnergyJ / releasedCombustionEnergyJ) * 100));
    else addLine('Released Chemical Energy:      N/A (Calculation Error or Zero)'); addLine('  -> Energy percentages cannot be calculated.'); end

    % === Save summary to TXT ===
    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    safeCaseName = matlab.lang.makeValidName(caseNameToDisplay);
    outputFolderTxt = fullfile('output', 'Simulation_Summaries');
    if ~exist(outputFolderTxt, 'dir'), try mkdir(outputFolderTxt); catch ME_mkdir, warning('Could not create output directory "%s": %s.', outputFolderTxt, ME_mkdir.message); outputFolderTxt = pwd; end; end
    summaryFilename = sprintf('%s_summary_%s.txt', safeCaseName, timestamp);
    summaryFullPath = fullfile(outputFolderTxt, summaryFilename);

    fid = fopen(summaryFullPath, 'w');
    if fid == -1, warning('Could not open file "%s" for writing.', summaryFullPath);
    else
        % Write Header Info
        % (Keep existing ASCII art write code)
        asciiArt = ["  _______              _____ _           _   "; "|__   __|            / ____| |         | |  "; "    | | ___  _ __ ___| (___ | |__   ___ | |_ "; "    | |/ _ \| '__/ __|\___ \| '_ \ / _ \| __|"; "    | | (_) | |  \__ \____) | | | | (_) | |_ "; "    |_|\___/|_|  |___/_____/|_| |_|\___/ \__|"; "                                             "; "         TORSHOT - Internal Ballistics Simulator"; "================================================="; "" ];
        for i = 1:length(asciiArt), fprintf(fid, '%s\n', asciiArt{i}); end
        fprintf(fid, 'Simulation Report for: %s\nGenerated: %s\n\n', caseNameToDisplay, datestr(now));

        % --- Try to Load and Write Component Data to TXT ---
        componentDataForTxt = []; componentDataLoadedForTxt = false;
        if canLoadComponentData && exist('loadComponentData', 'file') == 2
            fprintf(fid, '--- Input Component Data ---\n');
            try
                fprintf('Loading component data for TXT summary...\n');
                componentDataForTxt = loadComponentData(cartridgeName, bulletName, powderName);
                if isstruct(componentDataForTxt) && isfield(componentDataForTxt, 'cartridge')
                   componentDataLoadedForTxt = true;
                   printStructFields = @(fh, ds, t) printComponentStruct(fh, ds, t); % Local handle
                   printStructFields(fid, componentDataForTxt.cartridge, 'Cartridge Parameters');
                   printStructFields(fid, componentDataForTxt.bullet, 'Bullet Parameters');
                   printStructFields(fid, componentDataForTxt.powder, 'Powder Parameters');
                   % Also print loaded barrel data if available
                   if barrelDataLoaded
                       printStructFields(fid, barrelData, 'Barrel Parameters');
                   else
                       fprintf(fid, 'Barrel Parameters:\n  (Data not loaded or invalid)\n\n');
                   end
                else fprintf(fid, '  (Failed to load valid component data structure)\n'); end
            catch ME_load_txt, fprintf(fid, '  (Error loading component data: %s)\n', ME_load_txt.message); end
            fprintf(fid, '----------------------------\n\n');
        else fprintf(fid, '--- Input Component Data: NOT AVAILABLE (Cannot Load) ---\n\n'); end

        % --- Write Simulation Summary and Energy Balance ---
        fprintf(fid, '--- Simulation Results Summary ---\n'); % Add header
        for i = 1:length(summaryLines)
            if startsWith(summaryLines{i}, '---'), fprintf(fid, '%s\n', summaryLines{i}); else fprintf(fid, '  %s\n', summaryLines{i}); end
        end
        fclose(fid);
        fprintf('Summary saved to: %s\n', summaryFullPath);
    end

    % === Save .mat ===
    % (Component data loading logic - keep existing code)
    componentDataForMat = []; componentDataLoadedForMat = false;
    if canLoadComponentData && exist('loadComponentData', 'file') == 2
         try
             if componentDataLoadedForTxt && ~isempty(componentDataForTxt), componentDataForMat = componentDataForTxt; componentDataLoadedForMat = true; disp('Reusing component data loaded for TXT summary for MAT.');
             else componentDataForMat = loadComponentData(cartridgeName, bulletName, powderName); if isstruct(componentDataForMat) && isfield(componentDataForMat, 'cartridge'), componentDataLoadedForMat = true; disp('Successfully loaded component data for MAT file.'); else warning('Loaded component data for MAT file appears invalid.'); end; end
         catch ME_load_mat, warning('analyzeAndSaveResults: Failed to load componentData for MAT file: %s.', ME_load_mat.message); end
    end

    try
        outputFolder = fullfile('output', 'Matlab_Simulation_Data');
        if ~exist(outputFolder, 'dir'), try mkdir(outputFolder); catch ME_mkdir_mat, warning('Could not create .mat output directory "%s". Saving to pwd.', outputFolder); outputFolder = pwd; end; end
        matFilename = sprintf('%s_results_%s.mat', safeCaseName, timestamp);
        fullOutputPath = fullfile(outputFolder, matFilename);

        % --- Save all relevant data ---
        varsToSave = {'parameters', 'simulationResults'}; % Start with basic data
        % Add component data if loaded
        if componentDataLoadedForMat && ~isempty(componentDataForMat), varsToSave{end+1} = 'componentDataForMat'; end
        % Add barrel data if loaded
        if barrelDataLoaded, varsToSave{end+1} = 'barrelData'; end
        % Add stress results if calculated
        if ~isempty(stressResults), varsToSave{end+1} = 'stressResults'; end

        save(fullOutputPath, varsToSave{:}); % Save all specified variables
        fprintf('MAT data saved to: %s (Variables: %s)\n', fullOutputPath, strjoin(varsToSave, ', '));

    catch ME_savemat
        fprintf(2, 'Error saving .mat file: %s\n', ME_savemat.message);
    end

    % === Output struct ===
    % (Keep existing code)
    analysisOutput = struct();
    analysisOutput.finalVelocityMps = finalVelocitySim;
    analysisOutput.finalOmegaRpm = finalOmegaSimRpm;
    analysisOutput.maxPressurePa = maxPressureSimPa;
    analysisOutput.exitTimeS = exitTimeSimS;
    if initialPropMass > 1e-12, analysisOutput.propellantConsumedFraction = propellantConsumedKg / initialPropMass; else analysisOutput.propellantConsumedFraction = NaN; end
    analysisOutput.totalHeatLossJ = totalHeatLossJ;
    analysisOutput.totalFrictionWorkJ = totalFrictionWorkJ;
    % Add stress info to output struct
    if ~isempty(stressResults)
        analysisOutput.maxEquivalentStressPa = stressResults.max_eq_stress_Pa;
        analysisOutput.minSafetyFactor = stressResults.min_safety_factor;
    else
        analysisOutput.maxEquivalentStressPa = NaN;
        analysisOutput.minSafetyFactor = NaN;
    end


    disp('Analysis and saving completed.');

end % End of function analyzeAndSaveResults


% --- Helper function to print struct fields to file ---
% (Keep existing printComponentStruct function)
function printComponentStruct(fid, dataStruct, title)
    fprintf(fid, '%s:\n', title);
    if isempty(dataStruct) || ~isstruct(dataStruct), fprintf(fid, '  (Data not available)\n\n'); return; end
    fields = fieldnames(dataStruct); maxLen = 0;
    for k = 1:length(fields), maxLen = max(maxLen, length(fields{k})); end
    formatStr = sprintf('  %%-%ds : %%s\\n', maxLen + 2);
    for k = 1:length(fields)
        fn = fields{k}; v = dataStruct.(fn); vs = '';
        if ischar(v), vs = v;
        elseif isnumeric(v)
            if isscalar(v), if abs(v) > 1e5 || (abs(v) < 1e-4 && v ~= 0), vs = sprintf('%.4e', v); else vs = sprintf('%.4f', v); end
            else vs = sprintf('[%dx%d %s]', size(v,1), size(v,2), class(v)); end
        elseif islogical(v), if v, vs = 'true'; else vs = 'false'; end
        else vs = sprintf('(%s)', class(v)); end
        fprintf(fid, formatStr, fn, vs);
    end
    fprintf(fid, '\n');
end