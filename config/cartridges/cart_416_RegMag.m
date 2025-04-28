% =========================================================
% File: config/cartridges/RemMag_416.m
% Cartridge Data: .416 Remington Magnum
% Primary Units: Inches (in), Grains H2O (grH2O)
% Data Source: QuickLOAD V3.6 Screenshot (image_492922.png)
% Date Created: 2025-04-13
% =========================================================

disp('Executing config/cartridges/RemMag_416.m');

cartridgeData = struct();

% --- Identification ---
cartridgeData.cartridgeName = '.416 Remington Magnum (QL Data)';

% --- Dimensional and Capacity Data (from QuickLOAD Screenshot) ---
cartridgeData.caseLength_in        = 2.850;    % Case Length in INCHES
cartridgeData.maxCaseCapacity_grH2O = 107.00;   % Maximum Case Capacity (overflow) in GRAINS H2O
cartridgeData.boreDiameter_in      = 0.416;    % Bore Diameter (Groove Caliber) in INCHES

% --- Optional Data ---
cartridgeData.primerType = ''; % Example: 'Large Rifle Magnum' (Not specified in QL screenshot)

% --- Notes ---
cartridgeData.notes = 'Dimensional and capacity data sourced from QuickLOAD v3.6 screenshot (image_492922.png). Units: in, grH2O.';

disp('Struct "cartridgeData" defined for RemMag_416.');

% DO NOT add 'clear' at the end