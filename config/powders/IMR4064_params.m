% --- Parameter Definition for IMR 4064 ---
% NOTE: Incorporates data directly from QuickLOAD where possible.
%       Other parameters approximated from Varget data or based on assumptions.
%       Verify ALL approximations for your specific simulation needs.

imr4064PowderParams = struct(...
    'powderName', 'IMR 4064', ...
    'solidDensity_kg_m3', 1580, ... % Propellant Solid Density [kg/m^3] - FROM QuickLOAD
    'propellantDensity_rho_s', 1580, ... % Solid density used in ODEs [kg/m^3] - FROM QuickLOAD
    'impetus_F', 3880e3, ...        % Impetus (Heat of Explosion) [J/kg] - FROM QuickLOAD
    'covolume_b', 1e-04, ...        % Co-volume [m^3/kg] - APPROXIMATED (from Varget, QuickLOAD 'b' factor incompatible)
    'adiabaticFlameTemp_T_flame', 3400, ... % Adiabatic Flame Temperature [K] - APPROXIMATED (from Varget, not in QuickLOAD data)
    'specificHeatRatio_gamma', 1.2380, ... % Ratio Cp/Cv [-] - FROM QuickLOAD
    'molarMass_M_gas', 3.5e-02, ... % Average Molar Mass of gas products [kg/mol] - APPROXIMATED (from Varget, not in QuickLOAD data)
    'burnRateExponent_beta',   0.513488, ... % Combustion law exponent [-] - APPROXIMATED (from Varget, QuickLOAD uses different model: Ba/a0/z1)
    'burnRateCoeff_a', 7.259481e-06, ... % Combustion law coefficient [m/s / Pa^beta] - APPROXIMATED (from Varget, QuickLOAD uses different model)
    'specificSurfaceArea_sigma', 5.65, ... % Specific surface area [m^2/kg] - APPROXIMATED (from Varget, QuickLOAD uses different model)
    'formFunctionParam_theta', 3 ... % Form function parameter theta [-] - APPROXIMATED (Using QuickLOAD 'Factor b', assuming correlation observed with Varget holds)
); 

% --- Save to .mat File in 'powders' folder ---

% Desired subfolder name
subfolderName = fullfile('config', 'powders'); % Use fullfile for OS compatibility

% Filename to save
matFilename = 'IMR4064_params.mat';

% Check if the subfolder exists in the current directory
% If not, create it.
if ~exist(subfolderName, 'dir')
    fprintf('Folder "%s" not found in the current directory, creating it...\n', subfolderName);
    mkdir(subfolderName);
end

% Construct the full path to the file
fullPath = fullfile(subfolderName, matFilename);

% Save the struct to the .mat file specifying the full path
try
    % Save using the struct name as the variable name inside the MAT file
    save(fullPath, 'imr4064PowderParams');
    fprintf('File "%s" created successfully at path: "%s"\n', matFilename, fullPath);
catch ME
    fprintf('ERROR saving file to "%s":\n%s\n', fullPath, ME.message);
    fprintf('Check write permissions or the path.\n');
end
% Display the parameters (optional)
disp('Parameters for IMR 4064 (Combined QuickLOAD data and Varget approximations):');
disp(imr4064PowderParams);

% --- Code to save to .mat file (adapt as needed) ---
% ... (saving code similar to Varget_params.txt) ...