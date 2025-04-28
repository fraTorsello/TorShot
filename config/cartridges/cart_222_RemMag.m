% =========================================================
% File: config/cartridges/RemMag_222.m
% Cartridge Data: .222 Remington Magnum
% Primary Units: Inches (in), Grains H2O (grH2O)
% Data Source: QuickLOAD V3.6 Data
% Date Created: 2025-04-13
% =========================================================

disp('Executing config/cartridges/cart_222_RemMag.m');

cartridgeData = struct();

% --- Identification ---
cartridgeData.cartridgeName = 'cart_222_RemMag';

% --- Dimensional and Capacity Data ---
cartridgeData.caseLength_in = 1.850;    % Case Length in INCHES
cartridgeData.maxCaseCapacity_grH2O = 30.50;   % Maximum Case Capacity in GRAINS H2O
cartridgeData.boreDiameter_in = 0.224;    % Bore Diameter (Groove Caliber) in INCHES

% --- Optional Data ---
cartridgeData.primerType = ''; %  (Not specified in QL data)

% --- Notes ---
cartridgeData.notes = 'Dimensional and capacity data from QuickLOAD v3.6. Units: in, grH2O. Bore diameter is groove caliber.';

disp('Struct "cartridgeData" defined for RemMag_222.');

% No 'clear' command is needed at the end of a data definition script.
