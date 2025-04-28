% TORSHOT/gui/openPowderParamsCallback.m
function openPowderParamsCallback(hObject, eventdata)
% Apre la finestra parametri polvere, popolandola con i dati correnti.
% INPUT: hObject - Handle del controllo (pushbuttonPowderParams).
%        eventdata - Dati evento (non usati).
% NOTA: Questa funzione ora vive nel suo file .m e assume che logStatus
%       e guidata siano accessibili nel path.

    hFigure = ancestor(hObject, 'figure'); % Ottieni handle figura principale
    mainHandles = guidata(hFigure); % Ottieni handles figura principale

    % Controlla se logStatus esiste prima di chiamarla
    logFuncExists = exist('logStatus', 'file') == 2;
    if logFuncExists, logStatus('Edit Powder Params button pressed.'); else, disp('Edit Powder Params button pressed.'); end

    % Verifica se i componenti (e quindi la polvere) sono caricati
    if ~isfield(mainHandles,'isComponentsLoaded') || ~mainHandles.isComponentsLoaded || ...
       ~isfield(mainHandles,'componentData') || isempty(mainHandles.componentData) || ...
       ~isfield(mainHandles.componentData, 'powder') || isempty(mainHandles.componentData.powder)
        if logFuncExists, logStatus('Load components (including powder) first.'); else, disp('Load components (including powder) first.'); end
        uiwait(msgbox('Load components first to edit powder parameters.', 'Load Required', 'warn'));
        return;
    end

    % Verifica se l'handle della finestra polvere è valido
    if ~isfield(mainHandles, 'powderFigure') || isempty(mainHandles.powderFigure) || ~ishandle(mainHandles.powderFigure)
         msg = 'ERROR: Powder parameters window handle is invalid or missing.';
         if logFuncExists, logStatus(msg); else, disp(msg); end
         uiwait(msgbox('Error: Could not find the powder parameters window.', 'Window Error', 'error'));
         return;
    end

    try
        hPowderFig = mainHandles.powderFigure;
        powderHandles = guidata(hPowderFig); % Ottieni handles finestra polvere

        % Verifica se powderHandles è valido e contiene parameterFieldNames
        if isempty(powderHandles) || ~isfield(powderHandles, 'parameterFieldNames')
             msg = 'ERROR: Could not retrieve handles or field names from powder window.';
             if logFuncExists, logStatus(msg); else, disp(msg); end
             uiwait(msgbox('Error: Could not access powder window internals.', 'Window Error', 'error'));
             return;
        end

        loadedPowderData = mainHandles.componentData.powder; % Dati polvere caricati (SI)

        if logFuncExists, logStatus('Populating powder window with current parameters...'); else, disp('Populating powder window...'); end
        powderFieldNames = powderHandles.parameterFieldNames; % Nomi campi definiti in powderParametersGUI

        for i = 1:length(powderFieldNames)
             fieldName = powderFieldNames{i};
             editHandleName = ['edit', fieldName];

             if isfield(powderHandles, editHandleName) && ishandle(powderHandles.(editHandleName))
                 % Ottieni valore corrente dai dati caricati
                 currentValue = NaN;
                 if isfield(loadedPowderData, fieldName)
                     currentValue = loadedPowderData.(fieldName);
                 else
                      warnMsg = sprintf('WARN: Field "%s" not found in loaded powder data.', fieldName);
                      if logFuncExists, logStatus(warnMsg); else, disp(warnMsg); end
                 end

                 % Imposta il campo editabile nella finestra polvere
                 if ~isnan(currentValue)
                     % Formatta come stringa (usa esponenziale per numeri piccoli/grandi)
                     if abs(currentValue) > 1e5 || (abs(currentValue) < 1e-4 && currentValue ~= 0)
                          set(powderHandles.(editHandleName), 'String', sprintf('%.6e', currentValue));
                     else
                          set(powderHandles.(editHandleName), 'String', sprintf('%.6f', currentValue));
                     end
                 else
                     set(powderHandles.(editHandleName), 'String', ''); % Campo vuoto se NaN o mancante
                 end
                 set(powderHandles.(editHandleName), 'BackgroundColor', [1 1 1]); % Reset colore sfondo
             else
                  warnMsg = sprintf('WARN: Powder GUI handle %s missing.', editHandleName);
                  if logFuncExists, logStatus(warnMsg); else, disp(warnMsg); end
             end
        end

        % Rendi visibile la finestra polvere
        set(hPowderFig, 'Visible', 'on');
        figure(hPowderFig); % Portala in primo piano
        if logFuncExists, logStatus('Powder parameters window opened.'); else, disp('Powder parameters window opened.'); end

    catch ME_openPowder
        errMsg = ['ERROR opening or populating powder window: ', ME_openPowder.message];
        if logFuncExists, logStatus(errMsg); else, disp(errMsg); end
        fprintf(2, 'Powder Window Error Stack:\n'); disp(ME_openPowder.getReport);
        uiwait(msgbox(sprintf('Error accessing powder window:\n%s', ME_openPowder.message), 'Window Error', 'error'));
    end
end % Fine openPowderParamsCallback