% =========================================================
% File: config/bullets/Mil_FMJ_52_8gr_222_545x39.m
% Bullet Data: Mil FMJ 52.8gr (.222) - Likely 5.45x39mm type
% Primary Units: Grains (gr), Inches (in)
% Data Source: QuickLOAD V3.6 Screenshot (image_47b2b1.png)
% Date Created: 2025-04-13
% =========================================================
disp('Executing config/bullets/cal_222_52_8gr_Mil_FMJ_52.m');

bulletData = struct();

% --- Identification ---
bulletData.bulletName = 'cal_222_52_8gr_Mil_FMJ_52';

% --- Physical Data (from QuickLOAD Screenshot) ---
bulletData.mass_gr     = 52.8;     % Mass in GRAINS
bulletData.diameter_in = 0.222;    % Diameter in INCHES
bulletData.length_in   = 1.004;    % Length in INCHES

% --- Seating Depth (from QuickLOAD Screenshot) ---
% This is the length of the bullet seated inside the case mouth.
bulletData.seatingDepth_in = 0.574; % Overall Seating Depth in INCHES

% --- Ballistic Data (EXTERNAL SOURCE REQUIRED) ---
bulletData.bc_g1 = NaN; % BC G1 - REQUIRES EXTERNAL DATA
bulletData.formFactor_g1 = NaN; % Form Factor G1 - REQUIRES CALCULATION (i=SD/BC) or EXTERNAL DATA
% Note: QL doesn't show BC/Form Factor in this view.
% Sectional Density (SD) can be calculated: SD = (mass_gr / 7000) / diameter_in^2
% SD = (52.8 / 7000) / 0.222^2 = 0.153 approx.

% --- Notes ---
bulletData.notes = ['Dimensional/mass/seating depth data from QuickLOAD v3.6 screenshot (image_47b2b1.png). ',...
                    'BALLISTIC COEFFICIENT (BC) and FORM FACTOR (i) VALUES ARE MISSING and must be added from an external source. ', ...
                    'Bullet appears to be a 5.45x39mm type FMJ. ', ...
                    'Calculated SD approx 0.153.'];

disp('Struct "bulletData" defined for Mil_FMJ_52_8gr_222_545x39 (BC/Form Factor MISSING).');

% DO NOT add 'clear' at the end