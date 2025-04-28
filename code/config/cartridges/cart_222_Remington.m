% =========================================================
% File: config/cartridges/Rem_223.m
% Cartridge Data: .223 Remington
% Primary Units: Inches (in), Grains H2O (grH2O)
% Data Source: QuickLOAD V3.6 Data (image_be1b80.png) & SAAMI Specs
% Date Created: 2025-04-14
% =========================================================

disp('Executing config/cartridges/cart_222_Remington.m');

cartridgeData = struct();

% --- Identification ---
cartridgeData.cartridgeName = 'cart_222_Remington';

% --- Dimensional and Capacity Data ---
% Case length from SAAMI specification
cartridgeData.caseLength_in = 1.760;    % Case Length in INCHES (SAAMI Max)

% Max capacity from QuickLOAD data provided
cartridgeData.maxCaseCapacity_grH2O = 28.80;   % Maximum Case Capacity in GRAINS H2O

% Bore diameter based on common .223/5.56mm standard (groove diameter)
cartridgeData.boreDiameter_in = 0.224;    % Bore Diameter (Groove Caliber) in INCHES

% --- Optional Data ---
cartridgeData.primerType = 'Small Rifle'; % Typical primer for .223 Rem

% --- Notes ---
cartridgeData.notes = 'Capacity data from QuickLOAD screenshot (image_be1b80.png). Case length from SAAMI spec. Bore diameter is standard .224 groove.';

disp('Struct "cartridgeData" defined for Rem_223.');

% No 'clear' command is needed at the end of a data definition script.
