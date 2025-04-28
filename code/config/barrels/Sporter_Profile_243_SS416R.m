% config/barrels/Sporter_Profile_243_Win_SS416R.m
% =========================================================
% File Dati Canna: Profilo Sporter .243 Winchester in Acciaio Inox 416R (Stima)
% =========================================================
% Descrizione: Dati stimati per una canna di profilo sporter/caccia
%              in acciaio INOX 416R (o simile), adatta per il .243 Winchester.
%              Assume una lunghezza comune per fucili da caccia.
% Unità Base: SI (Metri, Pascal)
% Fonti Valori: Stime ingegneristiche, valori tipici per acciai inox da canna.
% Data Creazione: 2025-04-22
% Ultima Modifica: 2025-04-22 (Materiale cambiato a SS416R)
% =========================================================

disp('Esecuzione config/barrels/Sporter_Profile_243_Win_SS416R.m');

barrelData = struct();

% --- Identificazione ---
barrelData.barrelName = 'Sporter Profile .243 Win (Stainless 416R Estimate)';

% --- Geometria (SI Units) ---
% Diametro esterno vicino alla camera (stima per profilo sporter)
% Esempio: 27 mm = 0.027 m. Valore di ESEMPIO.
barrelData.outerDiameter_m = 0.028; % Diametro esterno [m]

% Diametro interno della canna (Groove Diameter per .243 Win -> 0.243 pollici)
% 0.243 in * 0.0254 m/in = 0.0061722 m.
barrelData.boreDiameter_m = 0.0061722; % Diametro alesaggio (interno) [m]

% Lunghezza totale della canna
% Esempio: 22 pollici (comune per caccia) = 22 * 0.0254 = 0.5588 m
barrelData.barrelLength_m = 0.5588; % Lunghezza canna [m]

% --- Materiale (Valori Stimati per Acciaio Inox 416R) ---
barrelData.material = 'Stainless Steel 416R (Estimated Properties)'; % MODIFICATO

% Limite di snervamento (Yield Strength) - Stima conservativa per 416R
% Esempio: 780 MPa = 780e6 Pa
barrelData.yieldStrength_Pa = 780e6; % Limite di snervamento [Pa] % MODIFICATO

% Modulo di Young (Elastic Modulus) - Valore tipico per acciaio inox
% Esempio: 200 GPa = 200e9 Pa
barrelData.youngsModulus_Pa = 200e9; % Modulo di Young [Pa] % MODIFICATO

% Coefficiente di Poisson - Valore tipico per acciaio inox
% Esempio: 0.28
barrelData.poissonsRatio = 0.28; % Coefficiente di Poisson [-] % MODIFICATO

% --- Dati per Calcolo Termico (Opzionale, INDICATIVI) ---
% Massa della Canna - Stima per profilo sporter 22" (densità simile a 4140)
% Esempio: ~1.1 kg? VALORE PURAMENTE INDICATIVO
barrelData.massKg = 1.1; % Massa [kg] (INDICATIVA)

% Calore Specifico - Valore tipico per acciaio inox
% Esempio: 470 J/(kg*K)
barrelData.specificHeatJkgK = 470; % Calore specifico [J/(kg*K)] (INDICATIVO) % MODIFICATO

% --- Note ---
barrelData.notes = ['Valori geometrici e di materiale stimati per una canna sporter .243 Win in INOX 416R. ',...
                    'Diametro esterno (0.027m) è una stima vicino alla camera. ',...
                    'Diametro interno (alesaggio/groove) impostato a 0.00617m (0.243"). ',...
                    'Lunghezza canna impostata a 0.5588m (22 pollici). ',...
                    'Yield Strength (780MPa), Modulo E (200GPa), Poisson (0.28) sono stime per SS 416R. ',... % MODIFICATO
                    'Massa e calore specifico sono INDICATIVI. ',...
                    'USARE DATI SPECIFICI DELL''ARMA REALE PER ANALISI ACCURATE.']; % MODIFICATO

disp('Struct "barrelData" definita per Sporter_Profile_243_Win_SS416R.');

% NON aggiungere 'clear' alla fine