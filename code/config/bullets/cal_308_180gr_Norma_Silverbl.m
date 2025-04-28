% =========================================================
% File Dati Proiettile: Norma Silverbl 180gr (.308) - SKU 67628
% Unità Primarie: Grani (gr), Pollici (in)
% Fonte Dati: Screenshot QuickLOAD V3.6, G1 BC fornito dall'utente
% =========================================================
disp('Esecuzione config/bullets/Norma_Silverbl_180gr_308.m');

bulletData = struct();

% --- Identificazione ---
bulletData.bulletName = '.308 Norma Silverbl 180gr (67628)';

% --- Dati Fisici (da QuickLOAD) ---
bulletData.mass_gr     = 180.0;    % Massa in GRANI
bulletData.diameter_in = 0.308;    % Diametro in POLLICI
bulletData.length_in   = 1.180;    % Lunghezza in POLLICI

bulletData.seatingDepth_in = 0.394; % Recommended seating depth in INCHES
% --- Dati Balistici (da sito shydasoutdoorcenter) ---
bulletData.bc_g1 = 0.615; % BC G1 (Aggiornato)
% Sectional Density (SD) = 0.230 (dal sito)
% Form Factor (i = SD/BC) = 0.230 / 0.365 = 0.6301...
bulletData.formFactor_g1 = 0.615; % Form Factor G1 (Calcolato/Aggiornato)

% --- Note Aggiuntive (da sito shydasoutdoorcenter) ---
% bulletData.sku = '16315'; % Già nel nome
% bulletData.qty_per_box = 50;
% bulletData.manufacturer = 'Nosler';
% bulletData.country = 'United States';
% bulletData.cannelure = false; % No
% bulletData.base_type = 'Flat';
% bulletData.profile = 'Spitzer';
% bulletData.type = 'Partition'; % Già nel nome
% bulletData.lead_free = false; % No
bulletData.notes = ['Dati dimensionali/massa da QL/Shydas. ', ...
                    'BC/SD da Shydas. Form Factor calcolato (i=SD/BC). ', ...
                    'Base: Flat, Profile: Spitzer, Cannelure: No, LeadFree: No.'];


disp('Struct "bulletData" definita per Norma_Silverbl_180gr_308.');

% NON aggiungere 'clear' alla fine