% =========================================================
% File Dati Proiettile: Nosler Partition SP 95gr (.243) - SKU 16315
% Unità Primarie: Grani (gr), Pollici (in)
% Fonte Dati: Screenshot QuickLOAD & shydasoutdoorcenter.com
% =========================================================
disp('Esecuzione config/bullets/Nosler_Partition_SP_95gr_243.m');

bulletData = struct();

% --- Identificazione ---
bulletData.bulletName = '.243 Nosler Partition SP 95gr (SKU 16315)'; % Aggiunto SKU

% --- Dati Fisici (Unità Comuni) ---
bulletData.mass_gr     = 95.0;     % Massa in GRANI (Confermato)
bulletData.diameter_in = 0.243;    % Diametro in POLLICI (Confermato)
bulletData.length_in   = 1.025;    % Lunghezza in POLLICI (Aggiornato da OAL sito)

% --- Dati Balistici (da sito shydasoutdoorcenter) ---
bulletData.bc_g1 = 0.365; % BC G1 (Aggiornato)
% Sectional Density (SD) = 0.230 (dal sito)
% Form Factor (i = SD/BC) = 0.230 / 0.365 = 0.6301...
bulletData.formFactor_g1 = 0.630; % Form Factor G1 (Calcolato/Aggiornato)

bulletData.seatingDepth_in =0.365; % Recommended seating depth in INCHES
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

disp('Struct "bulletData" definita per Nosler_Partition_SP_95gr_243 (aggiornata da sito).');

% NON aggiungere 'clear' alla fine