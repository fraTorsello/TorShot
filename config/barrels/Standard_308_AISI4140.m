% config/barrels/Standard_308_AISI4140.m
% =========================================================
% File Dati Canna: Profilo Standard .308 Win in AISI 4140
% =========================================================
% Descrizione: Dati rappresentativi per una canna di profilo
%              standard/caccia in acciaio AISI 4140 o simile,
%              adatta per calibri come il .308 Winchester.
% Unità Base: SI (Metri, Pascal)
% Fonti Valori: Stime ingegneristiche e valori tipici per acciai da canna.
% Data Creazione: 2025-04-22
% Ultima Modifica: 2025-04-22 (Aggiunto boreDiameter_m e barrelLength_m) % MODIFIED
% =========================================================

disp('Esecuzione config/barrels/Standard_308_AISI4140.m');

barrelData = struct();

% --- Identificazione ---
barrelData.barrelName = 'Standard Profile .308 Win (AISI 4140 Estimate)';

% --- Geometria (SI Units) ---
% Diametro esterno nella zona vicino alla camera (critica per Pmax)
% Esempio: 30 mm = 0.030 m. Questo è un valore di ESEMPIO e
% andrebbe misurato o ottenuto dalle specifiche dell'arma reale.
barrelData.outerDiameter_m = 0.030; % Diametro esterno [m]

% *** NUOVO CAMPO AGGIUNTO ***
% Diametro interno della canna (Groove Diameter per .308 Win)
% 7.82 mm = 0.00782 m. Valore tipico, verificare per canne specifiche.
barrelData.boreDiameter_m = 0.00782; % Diametro alesaggio (interno) [m]

% *** NUOVO CAMPO AGGIUNTO ***
% Lunghezza totale della canna
% Esempio: 24 pollici = 24 * 0.0254 = 0.6096 m
barrelData.barrelLength_m = 0.6096; % Lunghezza canna [m] % ADDED

% --- Materiale (Valori Tipici per Acciaio AISI 4140/4150 Bonificato) ---
barrelData.material = 'AISI 4140 CrMo Steel (Estimated Properties)';

% Limite di snervamento (Yield Strength) - Valore conservativo ma robusto
% Esempio: 850 MPa = 850e6 Pa
barrelData.yieldStrength_Pa = 850e6; % Limite di snervamento [Pa]

% Modulo di Young (Elastic Modulus) - Valore tipico per acciaio
% Esempio: 205 GPa = 205e9 Pa
barrelData.youngsModulus_Pa = 205e9; % Modulo di Young [Pa] (Opzionale, per deformazioni)

% Coefficiente di Poisson - Valore tipico per acciaio
% Esempio: 0.29
barrelData.poissonsRatio = 0.29; % Coefficiente di Poisson [-] (Opzionale, per deformazioni)

% --- Dati per Calcolo Termico (Opzionale, se si vuole DeltaT canna) ---
% Massa della Canna - Dipende fortemente da lunghezza e profilo!
% Esempio per canna da 24" profilo standard: ~1.5 kg? VALORE PURAMENTE INDICATIVO
barrelData.massKg = 1.5; % Massa [kg] (INDICATIVA)

% Calore Specifico - Valore tipico per acciaio
% Esempio: 480 J/(kg*K)
barrelData.specificHeatJkgK = 480; % Calore specifico [J/(kg*K)] (INDICATIVO)


% --- Note ---
% Note aggiornate per includere boreDiameter_m e barrelLength_m
barrelData.notes = ['Valori geometrici e di materiale stimati per una canna standard in AISI 4140/4150. ',...
                    'Diametro esterno (0.030m) è una stima vicino alla camera. ',...
                    'Diametro interno (alesaggio/groove) impostato a 0.00782m (tipico .308 Win). ',... % Nota Aggiornata
                    'Lunghezza canna impostata a 0.6096m (24 pollici). ',... % Nota Aggiornata
                    'Yield Strength (850MPa), Modulo E (205GPa), Poisson (0.29) sono valori tipici. ',...
                    'Massa e calore specifico sono INDICATIVI per calcolo DeltaT. ',...
                    'USARE DATI SPECIFICI DELL''ARMA REALE PER ANALISI ACCURATE.'];

disp('Struct "barrelData" definita per Standard_308_AISI4140 (con boreDiameter_m e barrelLength_m).');

% NON aggiungere 'clear' alla fine