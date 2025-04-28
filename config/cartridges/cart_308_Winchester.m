% =========================================================
% File Dati Bossolo/Cartuccia: .308 Winchester (da QuickLOAD)
% Unità Primarie: Pollici (in), Grani H2O (grH2O)
% Fonte Dati: Screenshot QuickLOAD V3.6
% =========================================================

disp('Esecuzione config/cartridges/cart_308_Winchester.m');

cartridgeData = struct();

% --- Identificazione ---
cartridgeData.cartridgeName = '.308 Winchester';

% --- Dati Dimensionali e Capacità (da QuickLOAD) ---
cartridgeData.caseLength_in        = 2.014;    % Lunghezza bossolo in POLLICI
cartridgeData.maxCaseCapacity_grH2O = 56.00;   % Capacità massima (overflow) in GRANI H2O
cartridgeData.boreDiameter_in      = 0.308;    % Diametro foratura canna (Groove Caliber) in POLLICI

% --- Dati Opzionali ---
cartridgeData.primerType = ''; % Esempio: 'Large Rifle' (Non fornito da QL screenshot)

% --- Note ---
cartridgeData.notes = 'Dati dimensionali e capacità da screenshot QuickLOAD v3.6. Unità: in, grH2O.';

disp('Struct "cartridgeData" definita per Winchester_308_QL.');

% NON aggiungere 'clear' alla fine