% config/barrels/Standard_Varmint_222_Rem.m
% =========================================================
% File Dati Canna: Profilo Standard/Varmint .222 Remington in AISI 4140
% =========================================================
% Descrizione: Dati stimati per una canna di profilo standard/varmint
%              in acciaio AISI 4140, adatta per il .222 Remington.
%              Assume una lunghezza comune per fucili varmint.
% Unità Base: SI (Metri, Pascal)
% Fonti Valori: Stime ingegneristiche, valori tipici per acciai e dimensioni cal. 222.
% Data Creazione: 2025-04-22
% Ultima Modifica: 2025-04-22
% =========================================================

disp('Esecuzione config/barrels/Standard_Varmint_222_Rem.m');

barrelData = struct();

% --- Identificazione ---
barrelData.barrelName = 'Standard/Varmint Profile .222 Rem (AISI 4140 Estimate)';

% --- Geometria (SI Units) ---
% Diametro esterno vicino alla camera (stima per profilo standard/varmint)
% Esempio: 28 mm = 0.028 m. Valore di ESEMPIO.
barrelData.outerDiameter_m = 0.028; % Diametro esterno [m]

% Diametro interno della canna (Groove Diameter per .222 Rem -> 0.224 pollici)
% 0.224 in * 0.0254 m/in = 0.0056896 m.
barrelData.boreDiameter_m = 0.0056896; % Diametro alesaggio (interno) [m]

% Lunghezza totale della canna
% Esempio: 24 pollici (comune per varmint) = 24 * 0.0254 = 0.6096 m
barrelData.barrelLength_m = 0.6096; % Lunghezza canna [m]

% --- Materiale (Valori Tipici per Acciaio AISI 4140 Bonificato) ---
barrelData.material = 'AISI 4140 CrMo Steel (Estimated Properties)';

% Limite di snervamento (Yield Strength) - Valore conservativo
% Esempio: 850 MPa = 850e6 Pa
barrelData.yieldStrength_Pa = 850e6; % Limite di snervamento [Pa]

% Modulo di Young (Elastic Modulus) - Valore tipico per acciaio
% Esempio: 205 GPa = 205e9 Pa
barrelData.youngsModulus_Pa = 205e9; % Modulo di Young [Pa]

% Coefficiente di Poisson - Valore tipico per acciaio
% Esempio: 0.29
barrelData.poissonsRatio = 0.29; % Coefficiente di Poisson [-]

% --- Dati per Calcolo Termico (Opzionale, INDICATIVI) ---
% Massa della Canna - Stima per profilo varmint 24"
% Esempio: ~1.3 kg? VALORE PURAMENTE INDICATIVO
barrelData.massKg = 1.3; % Massa [kg] (INDICATIVA)

% Calore Specifico - Valore tipico per acciaio
% Esempio: 480 J/(kg*K)
barrelData.specificHeatJkgK = 480; % Calore specifico [J/(kg*K)] (INDICATIVO)

% --- Note ---
barrelData.notes = ['Valori geometrici e di materiale stimati per una canna varmint/standard .222 Rem. ',...
                    'Diametro esterno (0.028m) è una stima vicino alla camera. ',...
                    'Diametro interno (alesaggio/groove) impostato a 0.00569m (0.224"). ',...
                    'Lunghezza canna impostata a 0.6096m (24 pollici). ',...
                    'Yield Strength (850MPa), Modulo E (205GPa), Poisson (0.29) sono valori tipici per AISI 4140. ',...
                    'Massa e calore specifico sono INDICATIVI. ',...
                    'USARE DATI SPECIFICI DELL''ARMA REALE PER ANALISI ACCURATE.'];

disp('Struct "barrelData" definita per Standard_Varmint_222_Rem.');

% NON aggiungere 'clear' alla fine