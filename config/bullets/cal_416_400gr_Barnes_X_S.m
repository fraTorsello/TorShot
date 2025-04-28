% =========================================================
% File: config/bullets/Barnes_X_S_400gr_416.m
% Bullet Data: Barnes-X Spitzer (SKU 41690) 400gr (.416)
% Primary Units: Grains (gr), Inches (in)
% Data Source: QuickLOAD V3.6 Screenshot (image_492239.png)
% Date Created: 2025-04-13
% =========================================================
disp('Executing config/bullets/Barnes_X_S_400gr_416.m');

bulletData = struct();

% --- Identification ---
% Note: "S" likely stands for Spitzer. SKU added.
bulletData.bulletName = '.416 Barnes X S 400gr (SKU 41690)';

% --- Physical Data (from QuickLOAD Screenshot) ---
bulletData.mass_gr     = 400.0;    % Mass in GRAINS
bulletData.diameter_in = 0.416;    % Diameter in INCHES
bulletData.length_in   = 1.595;    % Length in INCHES

% --- Seating Depth (from QuickLOAD Screenshot) ---
% This is the length of the bullet seated inside the case mouth.
bulletData.seatingDepth_in = 0.795; % Overall Seating Depth in INCHES

% --- Ballistic Data (EXTERNAL SOURCE REQUIRED) ---
bulletData.bc_g1 = NaN; % BC G1 - REQUIRES EXTERNAL DATA (e.g., Barnes Website/Manual)
bulletData.formFactor_g1 = NaN; % Form Factor G1 - REQUIRES CALCULATION (i=SD/BC) or EXTERNAL DATA
% Note: QL doesn't show BC/Form Factor in this view.
% Sectional Density (SD) can be calculated: SD = (mass_gr / 7000) / diameter_in^2
% SD = (400.0 / 7000) / 0.416^2 = 0.330 approx.

% ---