% config/barrels/Safari_Profile_416_RemMag.m
% =========================================================
% File Dati Canna: Profilo Safari/Heavy .416 Rem Mag in AISI 4140
% =========================================================
% Descrizione: Dati stimati per una canna di profilo pesante/safari
%              in acciaio AISI 4140, adatta per il .416 Remington Magnum.
%              Assume una lunghezza comune per fucili da caccia grossa.
% Unità Base: SI (Metri, Pascal)
% Fonti Valori: Stime ingegneristiche, valori tipici per acciai e dimensioni cal. 416.
% Data Creazione: 2025-04-22
% Ultima Modifica: 2025-04-22
% =========================================================

disp('Esecuzione config/barrels/Safari_Profile_416_RemMag.m');

barrelData = struct();

% --- Identificazione ---
barrelData.barrelName = 'Safari/Heavy Profile .416 RemMag (AISI 4140 Estimate)';

% --- Geometria (SI Units) ---
% Diametro esterno vicino alla camera (stima per profilo pesante)
% Esempio: 32 mm = 0.032 m. Valore di ESEMPIO.
barrelData.outerDiameter_m = 0.032; % Diametro esterno [m]

% Diametro interno della canna (Groove Diameter per .416 -> 0.416 pollici)
% 0.416 in * 0.0254 m/in = 0.0105664 m.
barrelData.boreDiameter_m = 0.0105664; % Diametro alesaggio (interno) [m]

% Lunghezza totale della canna
% Esempio: 24 pollici (comune per safari) = 24 * 0.0254 = 0.6096 m
barrelData.barrelLength_m = 0.6096; % Lunghezza canna [m]

% --- Materiale (Valori Tipici per Acciaio AISI 4140 Bonificato) ---
barrelData.material = 'AISI 4140 CrMo Steel (Estimated Properties)';

% Limite di snervamento (Yield Strength) - Valore conservativo
% Esempio: 850 MPa = 850e6 Pa (potrebbe essere maggiore per canne magnum)
barrelData.yieldStrength_Pa = 850e6; % Limite di snervamento [Pa]

% Modulo di Young (Elastic Modulus) - Valore tipico per acciaio
% Esempio: 205 GPa = 205e9 Pa
barrelData.youngsModulus_Pa = 205e9; % Modulo di Young [Pa]

% Coefficiente di Poisson - Valore tipico per acciaio
% Esempio: 0.29
barrelData.poissonsRatio = 0.29; % Coefficiente di Poisson [-]

% --- Dati per Calcolo Termico (Opzionale, INDICATIVI) ---
% Massa della Canna - Stima per profilo pesante 24"
% Esempio: ~1.8 kg? VALORE PURAMENTE INDICATIVO
barrelData.massKg = 1.8; % Massa [kg] (INDICATIVA)

% Calore Specifico - Valore tipico per acciaio
% Esempio: 480 J/(kg*K)
barrelData.specificHeatJkgK = 480; % Calore specifico [J/(kg*K)] (INDICATIVO)

% --- Note ---
barrelData.notes = ['Valori geometrici e di materiale stimati per una canna pesante/safari .416 Rem Mag. ',...
                    'Diametro esterno (0.032m) è una stima vicino alla camera. ',...
                    'Diametro interno (alesaggio/groove) impostato a 0.01057m (0.416"). ',...
                    'Lunghezza canna impostata a 0.6096m (24 pollici). ',...
                    'Yield Strength (850MPa), Modulo E (205GPa), Poisson (0.29) sono valori tipici per AISI 4140. ',...
                    'Massa e calore specifico sono INDICATIVI. ',...
                    'USARE DATI SPECIFICI DELL''ARMA REALE PER ANALISI ACCURATE.'];

disp('Struct "barrelData" definita per Safari_Profile_416_RemMag.');

% NON aggiungere 'clear' alla fine