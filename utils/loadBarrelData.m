% utils/loadBarrelData.m
% (Assumes getfield_safe.m is now in the 'utils' folder and on the path)

function barrelDataSI = loadBarrelData(barrelSelection)
%loadBarrelData Loads barrel data from a specified .m file in config/barrels
%   and returns the data in SI units. Validates required fields.
%   MODIFIED: Added loading and validation for barrelLength_m.
%
%   INPUT:
%     barrelSelection: String containing the base name of the barrel
%                      configuration file (e.g., 'Standard_308_AISI4140').
%   OUTPUT:
%     barrelDataSI: Struct containing validated barrel data in SI units.

    if nargin < 1 || ~ischar(barrelSelection) || isempty(barrelSelection)
        error('loadBarrelData: barrelSelection input must be a non-empty string.');
    end

    fprintf('--- Loading Barrel Data: %s ---\n', barrelSelection);

    % --- Define Path ---
    configDir = 'config';
    BARREL_FOLDER = fullfile(configDir, 'barrels');

    % --- Build File Path ---
    barrelFile = fullfile(BARREL_FOLDER, [barrelSelection, '.m']);
    barrelScriptName = [barrelSelection, '.m'];

    % --- Check File Existence ---
    if ~exist(barrelFile, 'file')
        altBarrelFile = fullfile(pwd, BARREL_FOLDER, barrelScriptName);
         if ~exist(altBarrelFile, 'file')
            error('loadBarrelData: Barrel configuration file not found at "%s" or relative path.', barrelFile);
         else
            barrelFile = altBarrelFile;
            fprintf('Info: Using relative path for barrel file: %s\n', barrelFile);
         end
    end

    % --- Load Data ---
    loadedBarrelData = [];
    originalDir = pwd;
    try
        [folderPath, ~, ~] = fileparts(barrelFile);
        if ~isempty(folderPath) && exist(folderPath,'dir')
             fprintf('Changing directory to: %s\n', folderPath);
             cd(folderPath);
             cleanupObj = onCleanup(@() cd(originalDir)); % Ensure directory change back
        else
             error('loadBarrelData: Could not determine folder path for barrel file.');
        end

        fprintf('Executing barrel config script: %s\n', barrelScriptName);
        clear barrelData;
        run(barrelScriptName); % MUST define 'barrelData' struct

        if ~exist('barrelData', 'var') || ~isstruct(barrelData)
            error('loadBarrelData: Script "%s" did not define the struct "barrelData".', barrelScriptName);
        end
        loadedBarrelData = barrelData;
        disp('Barrel data loaded from .m file.');

    catch ME_load
        fprintf(2, 'Error loading/running barrel file "%s":\n%s\n', barrelFile, ME_load.message);
        fprintf(2, 'Error occurred in %s at line %d\n', ME_load.stack(1).file, ME_load.stack(1).line);
        rethrow(ME_load);
    end

    % --- Process and Validate SI Data ---
    disp('Validating barrel data (assuming SI units)...');
    barrelDataSI = struct();

    % Use getfield_safe (assuming it's now on the path in 'utils')
    barrelDataSI.barrelName       = getfield_safe(loadedBarrelData, 'barrelName', 'N/A');
    barrelDataSI.outerDiameter_m  = getfield_safe(loadedBarrelData, 'outerDiameter_m', NaN);
    barrelDataSI.boreDiameter_m   = getfield_safe(loadedBarrelData, 'boreDiameter_m', NaN);
    barrelDataSI.barrelLength_m   = getfield_safe(loadedBarrelData, 'barrelLength_m', NaN); % ADDED
    barrelDataSI.materialName     = getfield_safe(loadedBarrelData, 'material', 'Unknown');
    barrelDataSI.yieldStrength_Pa = getfield_safe(loadedBarrelData, 'yieldStrength_Pa', NaN);
    barrelDataSI.youngsModulus_Pa = getfield_safe(loadedBarrelData, 'youngsModulus_Pa', NaN);
    barrelDataSI.poissonsRatio    = getfield_safe(loadedBarrelData, 'poissonsRatio', NaN);
    barrelDataSI.notes            = getfield_safe(loadedBarrelData, 'notes', '');
    barrelDataSI.massKg            = getfield_safe(loadedBarrelData, 'massKg', NaN);
    barrelDataSI.specificHeatJkgK = getfield_safe(loadedBarrelData, 'specificHeatJkgK', NaN);


    % --- Validation ---
    isValid = true;
    errorMessages = {};

    if isnan(barrelDataSI.outerDiameter_m) || ~isnumeric(barrelDataSI.outerDiameter_m) || barrelDataSI.outerDiameter_m <= 0
        isValid = false; errorMessages{end+1} = 'outerDiameter_m missing, non-numeric, or non-positive.';
    end
    if isnan(barrelDataSI.boreDiameter_m) || ~isnumeric(barrelDataSI.boreDiameter_m) || barrelDataSI.boreDiameter_m <= 0
        isValid = false; errorMessages{end+1} = 'boreDiameter_m missing, non-numeric, or non-positive. (MUST be defined in the barrel .m file).';
    end
    % ADDED Validation for barrelLength_m
    if isnan(barrelDataSI.barrelLength_m) || ~isnumeric(barrelDataSI.barrelLength_m) || barrelDataSI.barrelLength_m <= 0
        isValid = false; errorMessages{end+1} = 'barrelLength_m missing, non-numeric, or non-positive. (Should be defined in barrel .m file).';
    end
     if isnan(barrelDataSI.yieldStrength_Pa) || ~isnumeric(barrelDataSI.yieldStrength_Pa) || barrelDataSI.yieldStrength_Pa <= 0
        isValid = false; errorMessages{end+1} = 'yieldStrength_Pa missing, non-numeric, or non-positive.';
    end

    if isValid && barrelDataSI.boreDiameter_m >= barrelDataSI.outerDiameter_m
         isValid = false;
         errorMessages{end+1} = sprintf('boreDiameter_m (%.4f) cannot be >= outerDiameter_m (%.4f).', barrelDataSI.boreDiameter_m, barrelDataSI.outerDiameter_m);
     end

    % Optional field warnings
    if ~isnan(barrelDataSI.youngsModulus_Pa) && (~isnumeric(barrelDataSI.youngsModulus_Pa) || barrelDataSI.youngsModulus_Pa <= 0)
         warning('loadBarrelData: youngsModulus_Pa found but is non-numeric or non-positive in %s.', barrelSelection);
         barrelDataSI.youngsModulus_Pa = NaN;
    elseif isnan(barrelDataSI.youngsModulus_Pa)
         fprintf('Info: Optional field youngsModulus_Pa not found or NaN in %s.\n', barrelSelection);
    end
     if ~isnan(barrelDataSI.poissonsRatio) && (~isnumeric(barrelDataSI.poissonsRatio) || barrelDataSI.poissonsRatio < 0 || barrelDataSI.poissonsRatio >= 0.5)
         warning('loadBarrelData: poissonsRatio found but is non-numeric or outside typical range [0, 0.5) in %s.', barrelSelection);
         barrelDataSI.poissonsRatio = NaN;
    elseif isnan(barrelDataSI.poissonsRatio)
         fprintf('Info: Optional field poissonsRatio not found or NaN in %s.\n', barrelSelection);
    end
     if ~isnan(barrelDataSI.massKg) && (~isnumeric(barrelDataSI.massKg) || barrelDataSI.massKg <= 0)
         warning('loadBarrelData: massKg found but is non-numeric or non-positive in %s.', barrelSelection);
         barrelDataSI.massKg = NaN;
     end
     if ~isnan(barrelDataSI.specificHeatJkgK) && (~isnumeric(barrelDataSI.specificHeatJkgK) || barrelDataSI.specificHeatJkgK <= 0)
         warning('loadBarrelData: specificHeatJkgK found but is non-numeric or non-positive in %s.', barrelSelection);
         barrelDataSI.specificHeatJkgK = NaN;
     end


    if ~isValid
        error('loadBarrelData: Validation failed for "%s":\n - %s', barrelSelection, strjoin(errorMessages, '\n - '));
    end

    fprintf('Barrel data validated successfully (SI units assumed from file).\n');
    fprintf('--- Finished Loading Barrel Data: %s ---\n', barrelSelection);

end % End function loadBarrelData