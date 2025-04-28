% =========================================================
% File Dati Bossolo/Cartuccia: .243 Winchester
% Unità Primarie: Pollici (in), Grani H2O (grH2O)
% Fonte Dati: Screenshot QuickLOAD
% =========================================================

disp('Esecuzione config/cartridges/Winchester_243.m');

cartridgeData = struct();

% --- Identificazione ---
cartridgeData.cartridgeName = '.243 Winchester';

% --- Dati Dimensionali e Capacità (Unità Comuni) ---
cartridgeData.caseLength_in        = 2.044;    % Lunghezza bossolo in POLLICI
cartridgeData.maxCaseCapacity_grH2O = 54.0;    % Capacità massima in GRANI H2O
cartridgeData.boreDiameter_in      = 0.243;    % Diametro foratura canna (Groove) in POLLICI

% --- Dati Opzionali ---
cartridgeData.primerType = ''; % Esempio: 'Large Rifle'

% --- Note ---
cartridgeData.notes = 'Dati dimensionali e capacità da screenshot QuickLOAD v3.6. Unità: in, grH2O.';

disp('Struct "cartridgeData" definita per Winchester_243.');

% NON aggiungere 'clear' alla fine