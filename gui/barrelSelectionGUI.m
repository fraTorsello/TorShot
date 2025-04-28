% gui/barrelSelectionGUI.m
% MODIFIED: Added 3D axes and dynamic plot update on popup selection.
% CORRECTED: Removed 'WindowStyle','modal' to allow dropdown interaction.
function hBarrelFig = barrelSelectionGUI(hMainFigure)
% Creates a separate GUI window for selecting the barrel configuration file,
% including a 3D preview.

    fprintf('\n--- Entering barrelSelectionGUI (with 3D Preview) ---\n');
    hBarrelFig = []; % Initialize output to empty

    % --- Basic Input Validation ---
    if nargin < 1 || ~ishandle(hMainFigure)
        error('barrelSelectionGUI: Valid main figure handle (hMainFigure) is required.');
    end

    % --- Define Figure Properties (Increased Size) ---
    figWidth = 550; figHeight = 450; % Increased size for plot
    try
        mainFigUnitsOrig = get(hMainFigure, 'Units');
        set(hMainFigure, 'Units', 'pixels');
        mainFigPosPixels = get(hMainFigure, 'Position');
        set(hMainFigure, 'Units', mainFigUnitsOrig);

        figLeft = mainFigPosPixels(1) + mainFigPosPixels(3) + 20;
        figBottom = mainFigPosPixels(2) + mainFigPosPixels(4) - figHeight - 50;
        screenSize = get(0,'ScreenSize');
        if figLeft + figWidth > screenSize(3) - 20, figLeft = mainFigPosPixels(1) - figWidth - 20; end
        if figBottom < 50, figBottom = 50; end
        if figLeft < 20, figLeft = 20; end
        if figBottom + figHeight > screenSize(4) - 50, figBottom = screenSize(4) - figHeight - 50; end
        figPositionPixels = [figLeft, figBottom, figWidth, figHeight];

    catch ME_pos
         fprintf(2, 'ERROR calculating figure position: %s\n', ME_pos.message);
         uiwait(msgbox('Error positioning barrel selection window.','Window Error','error'));
         return;
    end

    % --- Create Figure ---
    try
        hBarrelFig = figure('Name', 'Select Barrel & Preview', ...
                            'Units', 'pixels', ...
                            'Position', figPositionPixels, ...
                            'NumberTitle', 'off', ...
                            'MenuBar', 'none', ...
                            'ToolBar', 'figure', ...
                            'Visible', 'off', ...
                            'Resize', 'on', ... 
                            'Tag', 'barrelSelectionFigure');

        if ~ishandle(hBarrelFig), error('Figure creation failed.'); end
        set(hBarrelFig, 'Units', 'normalized'); % Switch units after setting position
        set(hBarrelFig, 'CloseRequestFcn', {@barrelFigureCloseRequest, hBarrelFig});

    catch ME_fig
        fprintf(2, 'ERROR creating figure window: %s\n', ME_fig.message);
        uiwait(msgbox(sprintf('Failed to create barrel selection window:\n%s', ME_fig.message),'Window Error','error'));
        if ishandle(hBarrelFig), try delete(hBarrelFig); catch, end; end
        hBarrelFig = []; return;
    end

    % --- Setup Handles ---
    try
        handles = struct();
        handles.hMainFigure = hMainFigure; % Store main figure handle
        handles.loadedBarrelData = [];    % Initialize field to store loaded data

        % --- Create Panels for Layout ---
        topPanelHeight = 0.22; % Height for selection controls + buttons
        handles.topPanel = uipanel('Parent', hBarrelFig, 'Title', '', 'BorderType','none', ...
                                   'Units', 'normalized', 'Position', [0, 1-topPanelHeight, 1, topPanelHeight]);
        handles.plotPanel = uipanel('Parent', hBarrelFig, 'Title', 'Barrel Preview', ...
                                    'FontWeight', 'bold', 'Units', 'normalized', ...
                                    'Position', [0.05, 0.05, 0.9, 1-topPanelHeight-0.07]);

    catch ME_panel
        fprintf(2, 'ERROR creating panels: %s\n', ME_panel.message);
        if ishandle(hBarrelFig), delete(hBarrelFig); end
        hBarrelFig = []; return;
    end

    % --- Create Controls inside Top Panel ---
    try
        % --- Selection Dropdown ---
        selX = 0.05; selY = 0.65; selW = 0.9; selH = 0.25;
        uicontrol('Parent', handles.topPanel, 'Style', 'text', 'String', 'Select Barrel:', ...
                  'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                  'Position', [selX, selY + selH*0.8, 0.3, selH*0.5]); % Label above
        handles.popupmenuBarrel = uicontrol('Parent', handles.topPanel, 'Style', 'popupmenu', ...
                                             'String', {'(Scanning...)'}, 'Units', 'normalized', ...
                                             'Position', [selX, selY, selW, selH], ...
                                             'Tag', 'popupmenuBarrel_selector', ...
                                             'Callback', {@popupmenuBarrel_Callback, hBarrelFig}); % ADDED CALLBACK

        % --- Buttons ---
        btnY = 0.15; btnH = 0.25; btnW = 0.3; gap = 0.05;
        btnX_OK = 0.5 - btnW - gap/2;
        btnX_Cancel = 0.5 + gap/2;
        handles.pushbuttonOK = uicontrol('Parent', handles.topPanel, 'Style', 'pushbutton', ...
                                         'String', 'Select', 'FontWeight', 'bold', 'Units','normalized', ...
                                         'Position', [btnX_OK, btnY, btnW, btnH], ...
                                         'Callback', {@barrelSelectOkayCallback, hBarrelFig});
        handles.pushbuttonCancel = uicontrol('Parent', handles.topPanel, 'Style', 'pushbutton', ...
                                         'String', 'Cancel', 'Units','normalized', ...
                                         'Position', [btnX_Cancel, btnY, btnW, btnH], ...
                                         'Callback', {@barrelFigureCloseRequest, hBarrelFig});
    catch ME_controls
        fprintf(2, 'ERROR creating controls: %s\n', ME_controls.message);
        if ishandle(hBarrelFig), delete(hBarrelFig); end
        hBarrelFig = []; return;
    end

    % --- Create Axes in Plot Panel ---
    try
        handles.axesBarrel3D = axes('Parent', handles.plotPanel, 'Units', 'normalized', ...
                                    'Position', [0.1 0.1 0.8 0.8], 'Tag', 'axesBarrel3D');
        title(handles.axesBarrel3D, 'Select barrel to preview');
        axis(handles.axesBarrel3D, 'off'); % Initially off
    catch ME_axes
         fprintf(2, 'ERROR creating 3D axes: %s\n', ME_axes.message);
        if ishandle(hBarrelFig), delete(hBarrelFig); end
        hBarrelFig = []; return;
    end

    % --- Store Handles ---
    try
        guidata(hBarrelFig, handles);
    catch ME_guidata
         fprintf(2, 'ERROR storing guidata: %s\n', ME_guidata.message);
         if ishandle(hBarrelFig), delete(hBarrelFig); end
         hBarrelFig = []; return;
    end

    % --- Populate Dropdown ---
    % --- Populate Dropdown ---
    try
        % --- Improved Path Finding ---
        baseDir = '';
        try
            % Try finding the main GUI file to get the base path reliably
            mainGuiPath = which('ballisticSimulatorGUI.m');
            if ~isempty(mainGuiPath)
                baseDir = fileparts(mainGuiPath); % Get the directory containing ballisticSimulatorGUI.m
                % Assuming ballisticSimulatorGUI.m is directly inside TORSHOT, baseDir IS the TORSHOT folder.
                 fprintf('DEBUG (barrelSelectionGUI): Found baseDir via `which`: %s\n', baseDir);
            else
                 % Fallback using mfilename or pwd if `which` fails
                 try baseDir = fileparts(fileparts(mfilename('fullpath'))); catch, baseDir = pwd; end % Go up one level from 'gui'
                 if isempty(baseDir) || strcmp(baseDir,'.'), baseDir = pwd; end
                 fprintf('DEBUG (barrelSelectionGUI): Using fallback baseDir: %s\n', baseDir);
            end
        catch ME_findBase
             baseDir = pwd; % Ultimate fallback
             fprintf('DEBUG (barrelSelectionGUI): Error finding baseDir (%s). Using pwd: %s\n', ME_findBase.message, baseDir);
        end

        barrelConfigFolder = fullfile(baseDir, 'config', 'barrels');
        fprintf('DEBUG (barrelSelectionGUI): Attempting to use barrel config folder: "%s"\n', barrelConfigFolder);

        % Check if the determined folder actually exists
        if ~exist(barrelConfigFolder, 'dir')
             error('Determined barrel configuration folder does not exist: %s', barrelConfigFolder);
        end
        % --- End Improved Path Finding ---

        if exist('populateSingleDropdown', 'file') == 2
            mainHandlesForLog = guidata(hMainFigure); % Get main handles for logging context
            if isempty(mainHandlesForLog)
                 warning('barrelSelectionGUI: Could not get guidata from main figure for logging.');
                 % Optionally, pass handles of the current figure instead,
                 % although populateSingleDropdown expects hMainFigure field
                 % For now, proceed but logging within populateSingleDropdown might fail.
            end

            fprintf('DEBUG (barrelSelectionGUI): Calling populateSingleDropdown with Folder="%s", Handles Valid=%d\n', barrelConfigFolder, ~isempty(mainHandlesForLog));
            populateSingleDropdown(handles.popupmenuBarrel, barrelConfigFolder, {'*.m'}, 'Barrel', mainHandlesForLog); % Pass main GUI handles for logging
            fprintf('DEBUG (barrelSelectionGUI): populateSingleDropdown call finished.\n');
            set(handles.popupmenuBarrel, 'Value', 1); % Reset value
        else
            warning('barrelSelectionGUI: populateSingleDropdown function not found. Cannot populate list.');
            set(handles.popupmenuBarrel, 'String', {'Error: Func missing'});
        end
    catch ME_populate
         fprintf(2, 'ERROR during populateSingleDropdown call or path setup: %s\n', ME_populate.message);
         fprintf(2, 'Stack trace:\n'); disp(ME_populate.getReport);
         uiwait(msgbox(sprintf('Error populating barrel list:\n%s',ME_populate.message),'Popup Error','error'));
         if ishandle(hBarrelFig), set(handles.popupmenuBarrel, 'String', {'Error'}); end
    end


end % --- End of main barrelSelectionGUI function ---


% --- CALLBACKS and Helpers (Keep as before) ---

% --- NEW CALLBACK for Popup Menu ---
function popupmenuBarrel_Callback(hObject, ~, hBarrelFig)
    if ~ishandle(hBarrelFig), return; end
    barrelHandles = guidata(hBarrelFig);
    if isempty(barrelHandles), return; end

    try
        selectedIdx = get(hObject, 'Value');
        stringList = get(hObject, 'String');

        if ~iscell(stringList) || selectedIdx <= 1 % Index 1 is the placeholder
            cla(barrelHandles.axesBarrel3D);
            title(barrelHandles.axesBarrel3D, 'Select barrel to preview');
            axis(barrelHandles.axesBarrel3D, 'off');
            barrelHandles.loadedBarrelData = []; % Clear loaded data
            guidata(hBarrelFig, barrelHandles);
            return;
        end

        selectedBarrelName = stringList{selectedIdx};
        if contains(selectedBarrelName, {'(','No ','Err'}) % Check for invalid entries
             cla(barrelHandles.axesBarrel3D); title(barrelHandles.axesBarrel3D, 'Invalid selection'); axis(barrelHandles.axesBarrel3D, 'off');
             barrelHandles.loadedBarrelData = []; guidata(hBarrelFig, barrelHandles); return;
        end

        % Load the selected barrel data
        if exist('loadBarrelData', 'file') == 2
            tempBarrelData = loadBarrelData(selectedBarrelName);
            if isstruct(tempBarrelData) && ~isempty(fieldnames(tempBarrelData))
                barrelHandles.loadedBarrelData = tempBarrelData; % Store loaded data
                guidata(hBarrelFig, barrelHandles); % Save handles with new data

                % Plot the 3D barrel model
                if exist('plotBarrel3D', 'file') == 2
                    plotBarrel3D(barrelHandles.axesBarrel3D, tempBarrelData); % Call the 3D plot function
                else
                    title(barrelHandles.axesBarrel3D, 'plotBarrel3D not found'); axis(barrelHandles.axesBarrel3D, 'on'); grid(barrelHandles.axesBarrel3D, 'on');
                end
            else
                cla(barrelHandles.axesBarrel3D); title(barrelHandles.axesBarrel3D, 'Error loading data'); axis(barrelHandles.axesBarrel3D, 'off');
                barrelHandles.loadedBarrelData = []; guidata(hBarrelFig, barrelHandles);
            end
        else
             cla(barrelHandles.axesBarrel3D); title(barrelHandles.axesBarrel3D, 'loadBarrelData not found'); axis(barrelHandles.axesBarrel3D, 'off');
             barrelHandles.loadedBarrelData = []; guidata(hBarrelFig, barrelHandles);
        end

    catch ME_popup
        fprintf(2, 'ERROR in popupmenuBarrel_Callback: %s\n', ME_popup.message);
        try cla(barrelHandles.axesBarrel3D); title(barrelHandles.axesBarrel3D, 'Error in Callback'); axis(barrelHandles.axesBarrel3D, 'off'); catch; end
        barrelHandles.loadedBarrelData = []; guidata(hBarrelFig, barrelHandles);
    end
end

% Callback OK Button
function barrelSelectOkayCallback(~, ~, hBarrelFig)
    if ~ishandle(hBarrelFig), return; end
    barrelHandles = guidata(hBarrelFig);
    hMainFigure = barrelHandles.hMainFigure;
    if ~ishandle(hMainFigure), warning('Barrel Select OK: Main figure handle invalid.'); try delete(hBarrelFig); catch, end; return; end
    mainHandles = guidata(hMainFigure);

    barrelDataToPass = barrelHandles.loadedBarrelData;

    if isempty(barrelDataToPass) || ~isstruct(barrelDataToPass)
         popup = barrelHandles.popupmenuBarrel; list = get(popup,'String'); idx = get(popup,'Value');
         if iscell(list) && idx > 1 && ~contains(list{idx},{'(','No ','Err'})
             try
                 if exist('loadBarrelData','file')==2, barrelDataToPass = loadBarrelData(list{idx}); else error('loadBarrelData not found.'); end
             catch ME_final_load, uiwait(msgbox(sprintf('Could not load final selection "%s":\n%s', list{idx}, ME_final_load.message), 'Load Error', 'error')); return; end
         else uiwait(msgbox('Please select a valid barrel configuration file.', 'Selection Required', 'warn')); return; end
    end

    if ~isstruct(barrelDataToPass) || isempty(fieldnames(barrelDataToPass)), uiwait(msgbox('Loaded barrel data is invalid. Cannot apply selection.', 'Load Error', 'error')); return; end

    try
        mainHandles.barrelData = barrelDataToPass;
        list = get(barrelHandles.popupmenuBarrel,'String'); idx = get(barrelHandles.popupmenuBarrel,'Value');
        if iscell(list) && idx > 1, mainHandles.selectedBarrelName = list{idx}; else mainHandles.selectedBarrelName = ''; end
        guidata(hMainFigure, mainHandles); % Save to main GUI handles
        localLogStatus(mainHandles, sprintf('Barrel "%s" data selected and stored successfully.', mainHandles.selectedBarrelName));

        if exist('updateUIState', 'file') == 2, updateUIState(hMainFigure); else warning('updateUIState function not found.'); end

        if isfield(mainHandles, 'updateMainBarrelPlot') && isa(mainHandles.updateMainBarrelPlot, 'function_handle')
             try mainHandles.updateMainBarrelPlot(hMainFigure, barrelDataToPass); localLogStatus(mainHandles, 'Updated 3D barrel plot in main GUI.'); catch ME_updateMainPlot, localLogStatus(mainHandles, ['Error updating main GUI barrel plot: ' ME_updateMainPlot.message]); end
        end

        delete(hBarrelFig); % Close selection window on success

    catch ME_apply
         localLogStatus(mainHandles, ['ERROR applying barrel selection: ', ME_apply.message]); fprintf(2, 'Barrel Apply Error Stack:\n'); disp(ME_apply.getReport); uiwait(msgbox(sprintf('Error applying barrel selection:\n%s', ME_apply.message), 'Apply Error', 'error'));
    end
end

% Callback Cancel/Close Button
function barrelFigureCloseRequest(~, ~, hBarrelFig)
    if ishandle(hBarrelFig)
        try
            barrelHandles = guidata(hBarrelFig);
            if isfield(barrelHandles, 'hMainFigure') && ishandle(barrelHandles.hMainFigure)
                mainHandles = guidata(barrelHandles.hMainFigure);
                localLogStatus(mainHandles, 'Barrel selection cancelled or window closed.');
                if isfield(mainHandles, 'barrelSelectionFigure') && ishandle(mainHandles.barrelSelectionFigure) && mainHandles.barrelSelectionFigure == hBarrelFig
                    mainHandles.barrelSelectionFigure = []; guidata(barrelHandles.hMainFigure, mainHandles);
                end
            end
        catch ME_close, fprintf(2, 'Error during barrel window close request: %s\n', ME_close.message); end
        delete(hBarrelFig);
    end
end

% --- Local Log Status Helper ---
function localLogStatus(handles, message)
    logBoxHandle = []; hMainFigure = [];
    if isfield(handles, 'figure') && isfield(handles, 'logEditBox'), logBoxHandle = handles.logEditBox; hMainFigure = handles.figure;
    elseif isfield(handles, 'hMainFigure')
         hMainFigure = handles.hMainFigure;
         if ishandle(hMainFigure), mainHandles = guidata(hMainFigure); if isfield(mainHandles,'logEditBox'), logBoxHandle = mainHandles.logEditBox; end; end
    end
    if ~isempty(logBoxHandle) && ishandle(logBoxHandle)
        try
            hLog = logBoxHandle; hLogFig = ancestor(hLog, 'figure'); timestamp = datestr(now, 'HH:MM:SS');
            currentLog = get(hLog, 'String'); if isempty(currentLog), currentLog = {}; elseif ~iscell(currentLog), currentLog = {currentLog}; end
            if ischar(message), messageLines = strsplit(message, '\n'); elseif iscell(message), messageLines = message; else messageLines = {'Log Err: Invalid msg type'}; end
            newLogEntries = cell(sum(~cellfun('isempty', strtrim(messageLines))), 1); count = 0;
            for i = 1:length(messageLines), lineTxt = strtrim(messageLines{i}); if ~isempty(lineTxt), count=count+1; newLogEntries{count} = sprintf('[%s] %s', timestamp, lineTxt); end; end
            updatedLog = [currentLog; newLogEntries(1:count)]; maxLines = 200;
            if length(updatedLog) > maxLines, updatedLog = updatedLog(end-maxLines+1:end); end
            set(hLog, 'String', updatedLog); set(hLog, 'Value', length(updatedLog));
            if ishandle(hLogFig) && strcmp(get(hLogFig,'Visible'), 'off'), set(hLogFig,'Visible','on'); end
            drawnow limitrate;
        catch ME_log, fprintf(2, 'Internal localLogStatus error: %s\n', ME_log.message); disp(['LOG (Error): ', message]); end
    else, disp(['LOG (No Win/Handles): ', message]); end
end