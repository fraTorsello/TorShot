% gui/populateDisplayFields.m
function populateDisplayFields(hFigure)
% Aggiorna i campi di testo (read-only) nel pannello Dimensioni
% usando i dati caricati da handles.componentData (che è in SI).
% Esegue conversioni da SI base a unità di visualizzazione (gr, in, grH2O, in^2).
% INPUT: hFigure - Handle della figura principale della GUI.

    handles = guidata(hFigure);
    logStatus('Populating display fields with Imperial/Common units...');

    % Definisci i campi da pulire/popolare (usa i NUOVI NOMI TAG)
    fieldsToClear = {'textBulletName', 'textBulletMass_gr', 'textBulletLength_in', 'textBulletDiameter_in', ...
                     'textCartridgeName', 'textCaseLength_in', 'textMaxCaseCapacity_grH2O', 'textBoreArea_in2', ...
                     'textCartridgeOAL_in', 'textVolumeOccupiedBullet_grH2O', 'textUseableCaseCapacity_grH2O', 'textBulletTravel_in'};

    if isempty(handles.componentData) || ~isstruct(handles.componentData) ...
            || ~isfield(handles.componentData, 'bullet') || ~isfield(handles.componentData, 'cartridge')
        logStatus('No valid component data loaded to display.');
        for i = 1:length(fieldsToClear)
             if isfield(handles, fieldsToClear{i}) && ishandle(handles.(fieldsToClear{i})), set(handles.(fieldsToClear{i}), 'String', '---'); end
        end
        return;
    end

    cd = handles.componentData; % Dati componenti (in SI)
    C = handles.CONV;           % Costanti di conversione

    try
        % --- Dati Proiettile (Converti SI -> Display Units) ---
        set(handles.textBulletName, 'String', getfield_safe(cd.bullet, 'bulletName', 'N/A'));
        mass_gr = getfield_safe(cd.bullet, 'mass_kg', NaN) * C.GRAINS_PER_KG; % kg -> gr
        set(handles.textBulletMass_gr, 'String', iff(isnan(mass_gr), '---', sprintf('%.1f', mass_gr))); % Mostra 1 decimale per grani
        length_in = getfield_safe(cd.bullet, 'length_m', NaN) * C.INCHES_PER_METER; % m -> in
        set(handles.textBulletLength_in, 'String', iff(isnan(length_in), '---', sprintf('%.3f', length_in))); % Mostra 3 decimali per pollici
        diam_in = getfield_safe(cd.bullet, 'diameter_m', NaN) * C.INCHES_PER_METER; % m -> in
        set(handles.textBulletDiameter_in, 'String', iff(isnan(diam_in), '---', sprintf('%.3f', diam_in)));

        % --- Dati Bossolo (Converti SI -> Display Units) ---
        set(handles.textCartridgeName, 'String', getfield_safe(cd.cartridge, 'cartridgeName', 'N/A'));
        caseLen_in = getfield_safe(cd.cartridge, 'caseLength_m', NaN) * C.INCHES_PER_METER; % m -> in
        set(handles.textCaseLength_in, 'String', iff(isnan(caseLen_in), '---', sprintf('%.3f', caseLen_in)));
        maxCap_grH2O = getfield_safe(cd.cartridge, 'maxCaseCapacity_m3', NaN) * C.GRH2O_PER_M3; % m^3 -> grH2O
        set(handles.textMaxCaseCapacity_grH2O, 'String', iff(isnan(maxCap_grH2O), '---', sprintf('%.2f', maxCap_grH2O)));

        % Calcola e mostra Area Foratura in in^2
        boreDiam_m = getfield_safe(cd.cartridge, 'boreDiameter_m', NaN);
        if isnan(boreDiam_m)
            set(handles.textBoreArea_in2, 'String', '---');
        else
            bore_area_in2 = pi * (boreDiam_m * C.INCHES_PER_METER / 2)^2; % m -> in prima dell'area
            set(handles.textBoreArea_in2, 'String', sprintf('%.5f', bore_area_in2));
        end

        % --- Pulisci Campi Calcolati ---
        set(handles.textCartridgeOAL_in, 'String', '---');
        set(handles.textVolumeOccupiedBullet_grH2O, 'String', '---');
        set(handles.textUseableCaseCapacity_grH2O, 'String', '---');
        set(handles.textBulletTravel_in, 'String', '---');

        logStatus('Display fields updated.');
    catch ME_populate
         logStatus(['ERROR populating display fields: ', ME_populate.message]);
         for i = 1:length(fieldsToClear) % Pulisci tutti i campi in caso di errore
             if isfield(handles, fieldsToClear{i}) && ishandle(handles.(fieldsToClear{i})), set(handles.(fieldsToClear{i}), 'String', 'ERROR'); end
         end
    end
end

% --- Funzioni helper locali usate sopra ---
function val = getfield_safe(s, field, default)
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end
function out = iff(condition, trueVal, falseVal)
    if condition, out = trueVal; else out = falseVal; end
end