% =========================================================
% NUOVA FUNZIONE: loadComponentData
% Carica i dati separati per bossolo, proiettile e polvere.
% Esegue le conversioni in unità SI base (kg, m, m^3, Pa...).
% Sostituisce la vecchia logica di loadCaseParameters.
% =========================================================
function componentData = loadComponentData(cartridgeSelection, bulletSelection, powderSelection)
% INPUTS:
%   cartridgeSelection: Nome base del file bossolo (es. 'Winchester_243')
%   bulletSelection:    Nome base del file proiettile (es. 'Nosler_Partition_SP_95gr_243')
%   powderSelection:    Nome base del file polvere (es. 'Varget_params')
% OUTPUT:
%   componentData: Struct contenente i dati elaborati in SI.
%                  Fields: .cartridge (struct), .bullet (struct), .powder (struct)

fprintf('--- Caricamento Dati Componenti ---\n');
fprintf('Bossolo: %s | Palla: %s | Polvere: %s\n', cartridgeSelection, bulletSelection, powderSelection);

% --- Costanti di Conversione ---
KG_PER_GRAIN = 1 / 15432.35835;
METERS_PER_INCH = 0.0254;
M3_PER_GRAIN_H2O = KG_PER_GRAIN / 1000.0; % Densità acqua approx 1000 kg/m^3

% --- Percorsi Cartelle ---
configDir = 'config'; % Assumendo che sia nella cartella 'config' relativa alla CWD
CARTRIDGE_FOLDER = fullfile(configDir, 'cartridges');
BULLET_FOLDER = fullfile(configDir, 'bullets');
POWDER_FOLDER = fullfile(configDir, 'powders');

% --- Inizializza Output ---
componentData = struct('cartridge', [], 'bullet', [], 'powder', []);
loadedCartridgeData = [];
loadedBulletData = [];
loadedPowderData = [];

% --- Carica Dati Bossolo ---
cartridgeFile = fullfile(CARTRIDGE_FOLDER, [cartridgeSelection, '.m']);
if ~exist(cartridgeFile, 'file')
    error('loadComponentData: File bossolo non trovato: %s', cartridgeFile);
end
try
    clear cartridgeData; % Assicura pulizia da run precedenti
    run(cartridgeFile); % Esegue lo script .m che definisce 'cartridgeData'
    if ~exist('cartridgeData', 'var') || ~isstruct(cartridgeData)
        error('Lo script del bossolo "%s" non ha definito la struct "cartridgeData".', cartridgeFile);
    end
    loadedCartridgeData = cartridgeData; % Salva la struct caricata
    disp('Dati bossolo caricati.');
catch ME_cartridge
    fprintf(2, 'Errore durante caricamento/esecuzione file bossolo "%s":\n%s\n', cartridgeFile, ME_cartridge.message);
    rethrow(ME_cartridge);
end

% --- Carica Dati Proiettile ---
bulletFile = fullfile(BULLET_FOLDER, [bulletSelection, '.m']);
if ~exist(bulletFile, 'file')
    error('loadComponentData: File proiettile non trovato: %s', bulletFile);
end
try
    clear bulletData;
    run(bulletFile); % Esegue lo script .m che definisce 'bulletData'
    if ~exist('bulletData', 'var') || ~isstruct(bulletData)
        error('Lo script del proiettile "%s" non ha definito la struct "bulletData".', bulletFile);
    end
    loadedBulletData = bulletData;
    disp('Dati proiettile caricati.');
catch ME_bullet
    fprintf(2, 'Errore durante caricamento/esecuzione file proiettile "%s":\n%s\n', bulletFile, ME_bullet.message);
    rethrow(ME_bullet);
end

% --- Carica Dati Polvere ---
powderFile = fullfile(POWDER_FOLDER, [powderSelection, '.mat']);
if ~exist(powderFile, 'file')
    error('loadComponentData: File polvere non trovato: %s', powderFile);
end
try
    powderWorkspace = load(powderFile); % Carica il .mat in una struct temporanea
    powderVarNames = fieldnames(powderWorkspace);
    if numel(powderVarNames) ~= 1
        warning('Il file polvere "%s" contiene %d variabili. Uso la prima: "%s".', ...
                powderFile, numel(powderVarNames), powderVarNames{1});
    end
    loadedPowderData = powderWorkspace.(powderVarNames{1}); % Estrae la struct dati polvere
    if ~isstruct(loadedPowderData)
        error('La variabile "%s" caricata dal file polvere non è una struct.', powderVarNames{1});
    end
    disp('Dati polvere caricati.');
catch ME_powder
    fprintf(2, 'Errore durante caricamento file polvere "%s":\n%s\n', powderFile, ME_powder.message);
    rethrow(ME_powder);
end

% --- Elaborazione e Conversione in SI ---
disp('Conversione dati componenti in unità SI...');

% 1. Dati Bossolo (Conversione e Assegnazione)
siCartridge = struct();
siCartridge.cartridgeName = getfield_safe(loadedCartridgeData, 'cartridgeName', 'N/A');
siCartridge.caseLength_m = getfield_safe(loadedCartridgeData, 'caseLength_in', NaN) * METERS_PER_INCH;
siCartridge.maxCaseCapacity_m3 = getfield_safe(loadedCartridgeData, 'maxCaseCapacity_grH2O', NaN) * M3_PER_GRAIN_H2O;
siCartridge.boreDiameter_m = getfield_safe(loadedCartridgeData, 'boreDiameter_in', NaN) * METERS_PER_INCH;
siCartridge.primerType = getfield_safe(loadedCartridgeData, 'primerType', '');
siCartridge.notes = getfield_safe(loadedCartridgeData, 'notes', '');
% Validazione base
if isnan(siCartridge.caseLength_m) || siCartridge.caseLength_m <= 0 || ...
   isnan(siCartridge.maxCaseCapacity_m3) || siCartridge.maxCaseCapacity_m3 <= 0 || ...
   isnan(siCartridge.boreDiameter_m) || siCartridge.boreDiameter_m <= 0
    error('Dati bossolo caricati (%s) contengono valori mancanti, non validi o non positivi dopo la conversione.', cartridgeSelection);
end
componentData.cartridge = siCartridge;
fprintf('  Bossolo convertito in SI.\n');

% 2. Dati Proiettile (Conversione e Assegnazione)
siBullet = struct();
siBullet.bulletName = getfield_safe(loadedBulletData, 'bulletName', 'N/A');
siBullet.mass_kg = getfield_safe(loadedBulletData, 'mass_gr', NaN) * KG_PER_GRAIN;
siBullet.diameter_m = getfield_safe(loadedBulletData, 'diameter_in', NaN) * METERS_PER_INCH;
siBullet.length_m = getfield_safe(loadedBulletData, 'length_in', NaN) * METERS_PER_INCH;
siBullet.bc_g1 = getfield_safe(loadedBulletData, 'bc_g1', NaN); % Già numerico
siBullet.formFactor_g1 = getfield_safe(loadedBulletData, 'formFactor_g1', NaN); % Già numerico
siBullet.notes = getfield_safe(loadedBulletData, 'notes', '');

% Load and convert seating depth from bullet file
siBullet.seatingDepth_m = getfield_safe(loadedBulletData, 'seatingDepth_in', NaN) * METERS_PER_INCH;
if isnan(siBullet.seatingDepth_m) || siBullet.seatingDepth_m <= 0
     error('Bullet data file (%s) is missing a valid positive seatingDepth_in value.', bulletSelection);
end
% ---> END ADDED SECTION <---

% Validazione base
if isnan(siBullet.mass_kg) || siBullet.mass_kg <= 0 || ...
   isnan(siBullet.diameter_m) || siBullet.diameter_m <= 0 || ...
   isnan(siBullet.length_m) || siBullet.length_m <= 0
    error('Dati proiettile caricati (%s) contengono valori mancanti, non validi o non positivi dopo la conversione.', bulletSelection);
end
% Calcola Momento d'Inerzia Assiale (Stima come cilindro solido) - POTREBBE ESSERE NECESSARIO UN VALORE PIU' ACCURATO!
radius_m = siBullet.diameter_m / 2;
% Ip = 0.5 * m * r^2 per cilindro solido
siBullet.inertia_Ip_kgm2 = 0.5 * siBullet.mass_kg * radius_m^2; % Stima!
fprintf('  -> ATTENZIONE: Momento d''inerzia proiettile (inertia_Ip_kgm2) stimato come cilindro (%.3e kg*m^2).\n', siBullet.inertia_Ip_kgm2);
componentData.bullet = siBullet;
fprintf('  Proiettile convertito in SI.\n');

% 3. Dati Polvere (Assegnazione - Assumiamo già in SI nel file .mat)
% Verifica che i nomi dei campi nel file .mat corrispondano a quelli attesi dal simulatore
siPowder = struct();
siPowder.powderName = getfield_safe(loadedPowderData, 'powderName', powderSelection); % Usa nome file se non specificato
% Nomi campi come atteso da simulationOdes (o dalla struct parameters finale)
siPowder.propellantDensity_rho_s = getfield_safe(loadedPowderData, 'propellantDensity_rho_s', NaN); % kg/m^3
siPowder.impetus_F = getfield_safe(loadedPowderData, 'impetus_F', NaN); % J/kg
siPowder.covolume_b = getfield_safe(loadedPowderData, 'covolume_b', NaN); % m^3/kg
siPowder.adiabaticFlameTemp_T_flame = getfield_safe(loadedPowderData, 'adiabaticFlameTemp_T_flame', NaN); % K
siPowder.specificHeatRatio_gamma = getfield_safe(loadedPowderData, 'specificHeatRatio_gamma', NaN); % [-]
siPowder.molarMass_M_gas = getfield_safe(loadedPowderData, 'molarMass_M_gas', NaN); % kg/mol
siPowder.burnRateCoeff_a = getfield_safe(loadedPowderData, 'burnRateCoeff_a', NaN); % m/s / Pa^beta
siPowder.burnRateExponent_beta = getfield_safe(loadedPowderData, 'burnRateExponent_beta', NaN); % [-]
siPowder.specificSurfaceArea_sigma = getfield_safe(loadedPowderData, 'specificSurfaceArea_sigma', NaN); % m^2/kg
siPowder.formFunctionParam_theta = getfield_safe(loadedPowderData, 'formFunctionParam_theta', NaN); % [-]

% Calcola Costante Specifica Gas (se necessario)
if ~isnan(siPowder.molarMass_M_gas) && siPowder.molarMass_M_gas > 0
    R_univ = 8.314462;
    siPowder.specificGasConstant_R = R_univ / siPowder.molarMass_M_gas;
else
    siPowder.specificGasConstant_R = NaN;
     warning('Massa molare gas non valida, impossibile calcolare R_gas.');
end
% Validazione base (verifica solo alcuni campi chiave)
required_powder_fields = {'propellantDensity_rho_s', 'impetus_F', 'covolume_b', 'specificHeatRatio_gamma', 'burnRateCoeff_a', 'burnRateExponent_beta', 'specificSurfaceArea_sigma', 'formFunctionParam_theta', 'specificGasConstant_R'};
missing_powder = {};
for k=1:length(required_powder_fields)
    fld = required_powder_fields{k};
    if ~isfield(siPowder, fld) || isnan(siPowder.(fld))
        missing_powder{end+1} = fld;
    end
end
if ~isempty(missing_powder)
    error('Dati polvere (%s) mancano di campi essenziali o contengono NaN: %s', powderSelection, strjoin(missing_powder, ', '));
end
componentData.powder = siPowder;
fprintf('  Dati polvere verificati (assunti già in SI).\n');

fprintf('--- Caricamento e conversione SI completati ---\n');

end % Fine funzione loadComponentData

% --- Helper Function Interna ---
function val = getfield_safe(s, field, default)
    % Restituisce il valore del campo o un default se non esiste o è vuoto.
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end