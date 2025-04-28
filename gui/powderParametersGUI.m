function hPowderFig = powderParametersGUI(hMainFigure)
% Creates the GUI window for editing powder parameters.
% INPUT: hMainFigure - Handle to the main ballisticSimulatorGUI figure.
% OUTPUT: hPowderFig - Handle to the created (but initially invisible) powder figure.

fprintf('Creating Powder Parameters GUI...\n');

% --- Define Powder Field Names ---
powderFieldNames = { ...
    'propellantDensity_rho_s', 'covolume_b', 'adiabaticFlameTemp_T_flame', ...
    'specificHeatRatio_gamma', 'molarMass_M_gas', 'burnRateCoeff_a', ...
    'burnRateExponent_beta', 'formFunctionParam_theta', ...
    'specificSurfaceArea_sigma', 'impetus_F' ...
    };
powderHandles.parameterFieldNames = powderFieldNames;

% --- Create Figure ---
powderPosition = [0.7 0.5 0.25 0.45];
hPowderFig = figure('Name', 'Powder Parameters', ...
    'Units', 'normalized', 'Position', powderPosition, ...
    'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'Resize', 'on', 'Visible', 'off', ...
    'CloseRequestFcn', @powderFigureCloseRequest, ...
    'Tag', 'powderFigure');

% --- Create Panel ---
powderHandles.powderPanel = uipanel('Parent', hPowderFig, 'Title', 'Propellant Parameters (Edit Cautiously!)', ...
    'FontWeight', 'bold', 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);

% --- Create Controls ---
leftMargin = 0.05; labelWidth = 0.55; editWidth = 0.35; vSpacing = 0.01;

% Define Margins and Button Area Height (Normalized within the panel)
topMargin = 0.05; % Space above the first field
buttonHeight = 0.06;
spaceBelowButton = 0.01;
spaceAboveButton = 0.02; % << ADJUST THIS VALUE FOR MORE/LESS SPACE >>
buttonAreaTotalHeight = buttonHeight + spaceBelowButton + spaceAboveButton; % Total space reserved at the bottom

% Calculate available height for fields
topYBoundary = 1.0 - topMargin;
bottomYBoundary = buttonAreaTotalHeight;
availableHeightForFields = topYBoundary - bottomYBoundary;

% Calculate height per field, ensuring minimum visibility
numFields = length(powderFieldNames);
totalSpacingBetweenFields = max(0, numFields - 1) * vSpacing; % Spacing only exists between fields
heightPurelyForFields = availableHeightForFields - totalSpacingBetweenFields;

if numFields > 0 && heightPurelyForFields > 0
    fieldHeight = max(0.04, heightPurelyForFields / numFields); % Ensure min height
else
    fieldHeight = 0.05; % Default height if no fields or no space
end

currentY = topYBoundary; % Start placing from the top margin


% --- Nested Helper function ---
    function addPowderParamEdit(parentPanel, fieldName, labelText)
        % Label
        powderHandles.(['label', fieldName]) = uicontrol('Parent', parentPanel, 'Style', 'text', ...
            'String', [labelText, ':'], 'HorizontalAlignment', 'right', ...
            'Units', 'normalized', 'Position', [leftMargin, currentY - fieldHeight, labelWidth, fieldHeight]);
        % Edit Box
        powderHandles.(['edit', fieldName]) = uicontrol('Parent', parentPanel, 'Style', 'edit', ...
            'String', '', 'HorizontalAlignment', 'left', 'Tag', fieldName, 'Units', 'normalized', ...
            'Position', [leftMargin+labelWidth+0.01, currentY - fieldHeight, editWidth, fieldHeight], ...
            'Enable', 'on', 'BackgroundColor', [1 1 1]);
         currentY = currentY - fieldHeight - vSpacing;
    end % --- End of nested function ---

% Create all powder parameter fields
for k = 1:length(powderFieldNames)
    fieldName = powderFieldNames{k};
    labelText = regexprep(fieldName, '_', ' '); labelText = regexprep(labelText, '(\<[a-z])', '${upper($1)}');
    labelText = strrep(labelText,' rho s', ' (rho s)'); labelText = strrep(labelText, ' b', ' (b)'); labelText = strrep(labelText, ' T flame', ' (T flame)');
    labelText = strrep(labelText, ' gamma', ' (gamma)'); labelText = strrep(labelText, ' M gas', ' (M gas)'); labelText = strrep(labelText, ' a', ' (a)');
    labelText = strrep(labelText, ' beta', ' (beta)'); labelText = strrep(labelText, ' theta', ' (theta)');
    labelText = strrep(labelText, ' sigma', ' (sigma)'); labelText = strrep(labelText, ' F', ' (F)');
    addPowderParamEdit(powderHandles.powderPanel, fieldName, labelText);
end

% OK Button --- Positioned based on reserved space ---
buttonBottomY = spaceBelowButton; % Position relative to bottom margin
powderHandles.pushbuttonOK = uicontrol('Parent', powderHandles.powderPanel, 'Style', 'pushbutton', 'String', 'OK & Apply', 'FontWeight', 'bold', ...
    'Units','normalized', 'Position', [0.35, buttonBottomY, 0.3, buttonHeight], 'Callback', @powderOkayCallback);

% --- Store Handles ---
powderHandles.hMainFigure = hMainFigure;
guidata(hPowderFig, powderHandles);
fprintf('Powder Parameters GUI created.\n');

end % --- End of main powderParametersGUI function ---


%% --- CALLBACKS for Powder Window ---

function powderOkayCallback(hObject, ~)
hPowderFig = ancestor(hObject, 'figure'); powderHandles = guidata(hPowderFig); hMainFigure = powderHandles.hMainFigure;
if ~ishandle(hMainFigure), warning('Powder OK: Main figure handle invalid.'); set(hPowderFig, 'Visible', 'off'); return; end
mainHandles = guidata(hMainFigure); logStatus(hMainFigure, 'Applying powder parameters...');
tempParams = mainHandles.parametri; allValid = true; paramsModified = false; powderFieldNames = powderHandles.parameterFieldNames;
for i = 1:length(powderFieldNames)
     fieldName = powderFieldNames{i}; editHandleName = ['edit', fieldName];
     if isfield(powderHandles, editHandleName) && ishandle(powderHandles.(editHandleName))
        strValue = get(powderHandles.(editHandleName), 'String'); numValue = str2double(strValue);
        if isnan(numValue) || ~isreal(numValue), logStatus(hMainFigure, sprintf('POWDER ERR: Invalid number for "%s": "%s"', fieldName, strValue)); set(powderHandles.(editHandleName), 'BackgroundColor', [1 0.8 0.8]); allValid = false; continue; end
        validationPassed = true;
         if contains(fieldName,'Density') && numValue <= 0, validationPassed = false; logStatus(hMainFigure, sprintf('POWDER ERR: %s > 0',fieldName)); end
         if strcmp(fieldName,'specificHeatRatio_gamma') && numValue <= 1, validationPassed = false; logStatus(hMainFigure, sprintf('POWDER ERR: gamma > 1',fieldName)); end
         if strcmp(fieldName,'impetus_F') && numValue <= 0, validationPassed = false; logStatus(hMainFigure, sprintf('POWDER ERR: impetus_F > 0',fieldName)); end
        if ~validationPassed, set(powderHandles.(editHandleName), 'BackgroundColor', [1 0.8 0.8]); allValid = false; else
            currentValue = getParamValue(mainHandles.parametri, fieldName, NaN); if ~isequaln(currentValue, numValue), tempParams.(fieldName) = numValue; paramsModified = true; end; set(powderHandles.(editHandleName), 'BackgroundColor', [1 1 1]);
        end
     else logStatus(hMainFigure, sprintf('WARN: Powder GUI handle %s missing.', editHandleName)); allValid = false; end
end
if allValid
    if paramsModified, mainHandles.parametri = tempParams; guidata(hMainFigure, mainHandles); logStatus(hMainFigure, 'Powder params updated.'); checkFiConsistency(hMainFigure); else logStatus(hMainFigure, 'No powder values changed.'); end
    set(hPowderFig, 'Visible', 'off');
else logStatus(hMainFigure, 'Powder params NOT applied: errors.'); uiwait(msgbox('Invalid powder values entered.', 'Input Error', 'error')); end
end

function powderFigureCloseRequest(hObject, ~)
set(hObject, 'Visible', 'off');
try powderHandles = guidata(hObject); if isfield(powderHandles, 'hMainFigure') && ishandle(powderHandles.hMainFigure), logStatus(powderHandles.hMainFigure, 'Powder window closed.'); end; catch, end
end

function value = getParamValue(paramStruct, fieldName, defaultValue)
    if isfield(paramStruct, fieldName), value = paramStruct.(fieldName); if isempty(value), value = defaultValue; end; else value = defaultValue; end
end

function logStatus(hFigure, message)
    if ~ishandle(hFigure), fprintf('Log Err: Invalid fig handle.\n'); disp(message); return; end; try handles = guidata(hFigure); catch ME_guidata, fprintf('Log Err: guidata failed: %s\n', ME_guidata.message); disp(message); return; end; if ~isfield(handles, 'editStatusLog') || ~ishandle(handles.editStatusLog), fprintf('Log Err: Status log handle invalid.\n'); disp(message); return; end
    try timestamp = datestr(now, 'HH:MM:SS'); currentLog = get(handles.editStatusLog, 'String'); if isempty(currentLog), currentLog = {}; elseif ~iscell(currentLog), currentLog = {currentLog}; end; if ischar(message), messageLines = strsplit(message, '\n'); elseif iscell(message), messageLines = message; else messageLines = {'Log Err: Invalid msg type'}; end; newLogEntries = {}; for i = 1:length(messageLines), lineTxt = strtrim(messageLines{i}); if ~isempty(lineTxt), newLogEntries{end+1,1} = sprintf('[%s] %s', timestamp, lineTxt); end; end; updatedLog = [currentLog; newLogEntries]; maxLines = 150; if length(updatedLog) > maxLines, updatedLog = updatedLog(end-maxLines+1:end); end; set(handles.editStatusLog, 'String', updatedLog); set(handles.editStatusLog, 'Value', length(updatedLog)); drawnow limitrate;
    catch ME_log, fprintf(2, 'Internal logStatus error: %s\n', ME_log.message); fprintf(2, 'Original Log: '); if ischar(message), fprintf(2, '%s\n', message); elseif iscell(message), disp(message); end; end
end

function checkFiConsistency(hFigure)
    handles = guidata(hFigure); if ~isfield(handles,'parametri') || isempty(handles.parametri), return; end; params = handles.parametri; req_fields = {'impetus_F', 'adiabaticFlameTemp_T_flame', 'specificHeatRatio_gamma', 'specificGasConstant_R'}; pCheck = struct(); valid = true;
    for k = 1:length(req_fields), fn = req_fields{k}; if isfield(params, fn), v = params.(fn); else valid = false; break; end; if isnumeric(v) && isscalar(v) && ~isnan(v) && isreal(v) && v > 0, if strcmp(fn,'specificHeatRatio_gamma') && v <= 1, valid = false; break; end; pCheck.(fn) = v; else valid = false; break; end; end
    if valid, F_i=pCheck.impetus_F; T_0i=pCheck.adiabaticFlameTemp_T_flame; gamma_i=pCheck.specificHeatRatio_gamma; R_gas=pCheck.specificGasConstant_R; try Cv=R_gas/(gamma_i-1); Fi_chk=Cv*T_0i; rel_diff=abs(Fi_chk-F_i)/abs(F_i); tol=0.01; if rel_diff>tol, logStatus(hFigure,sprintf('WARN: F_i (%.3e) differs >%.0f%% from T_flame/gamma calc (%.3e).',F_i,tol*100,Fi_chk)); end; catch ME_check, logStatus(hFigure,sprintf('Info: Err during F_i consistency check: %s', ME_check.message)); end
    else logStatus(hFigure, 'Info: Cannot check F_i consistency (params missing/invalid).'); end
end