% utils/populateSingleDropdown.m
% Populates a single popup menu by scanning a folder.
% REVISED: Accepts handles from caller, uses hMainFigure for logging.
% CORRECTED: Updated debug print and logStatusHelper to check for 'figure' field.

function populateDropdown(hPopup, folder, extensions, typeName, handles_caller)
% Populates a single popup menu by scanning a folder.
% INPUTS:
%   hPopup:         Handle of the popup menu uicontrol.
%   folder:         Relative or absolute path to the folder to scan.
%   extensions:     Cell array of file extensions (e.g., {'*.m'}).
%   typeName:       String describing the type of item (e.g., 'Barrel', 'Powder').
%   handles_caller: Handles struct from the calling GUI (MUST contain the 'figure' field).

    % --- >>> CORRECTED DEBUG PRINT <<< ---
    fprintf('DEBUG (populateSingleDropdown %s): Received Folder="%s", hPopup Valid=%d, handles_caller Valid=%d\n', ...
            typeName, folder, ishandle(hPopup), isstruct(handles_caller));
    if isstruct(handles_caller) && isfield(handles_caller,'figure') % Check for 'figure' field
        fprintf('DEBUG (populateSingleDropdown %s): handles_caller.figure Valid=%d\n', typeName, ishandle(handles_caller.figure));
    elseif isstruct(handles_caller)
         fprintf('DEBUG (populateSingleDropdown %s): handles_caller struct DOES NOT contain ''figure'' field.\n', typeName); % Updated message
    else
         fprintf('DEBUG (populateSingleDropdown %s): handles_caller is NOT a struct.\n', typeName);
    end
    % --- >>> END CORRECTED DEBUG PRINT <<< ---


    % --- Input Validation (Corrected to check for 'figure') ---
    if nargin < 5 || ~ishandle(hPopup) || ~strcmp(get(hPopup,'Style'),'popupmenu') || ...
       ~ischar(folder) || isempty(folder) || ~iscell(extensions) || isempty(extensions) || ...
       ~ischar(typeName) || ~isstruct(handles_caller) || ~isfield(handles_caller, 'figure') || ~ishandle(handles_caller.figure) % Check for 'figure'
        warning('populateSingleDropdown: Invalid inputs (check handles_caller structure and figure field).'); % Updated warning
        try
            set(hPopup, 'String', {'Input Error'}, 'Value', 1, 'Enable', 'off');
        catch
            fprintf(2,'populateSingleDropdown: Failed to set Input Error on popup.\n');
        end
        return;
    end

    set(hPopup, 'String', {'Scanning...'}, 'Value', 1, 'Enable', 'off'); drawnow;

    fileList = {};
    foundFolder = false;
    folderPath = ''; % Initialize folderPath

    % --- Resolve Folder Path ---
    if exist(folder, 'dir')
        foundFolder = true;
        folderPath = folder;
    else
        % Use the logStatus helper defined below, passing handles_caller
        logStatusHelper(handles_caller, sprintf('ERROR: %s folder not found at: %s.', typeName, folder));
        set(hPopup, 'String', {['Folder Error']}, 'Enable', 'off', 'Value', 1);
        return;
    end

    % --- Scan Folder ---
    if foundFolder
        for k = 1:length(extensions)
            currentExtension = extensions{k};
            if ~startsWith(currentExtension, '*.')
                warning('populateDropdown: Extension format should be ''*.ext'', fixing: %s', currentExtension);
                currentExtension = ['*', currentExtension];
                if ~contains(currentExtension, '.'), currentExtension = [currentExtension, '.*']; end
                if ~startsWith(currentExtension, '*.'), currentExtension = ['*', currentExtension]; end
            end

            try
                files = dir(fullfile(folderPath, currentExtension));
                for i = 1:length(files)
                    if files(i).isdir, continue; end % Skip directories
                    [~, name, ~] = fileparts(files(i).name);
                    % Avoid hidden/system/temp files
                    if ~isempty(name) && ~startsWith(name, '.') && ~endsWith(name, '~')
                         fileList{end+1} = name; %#ok<AGROW>
                    end
                end
            catch ME_dir
                logStatusHelper(handles_caller, sprintf('Warning: Could not read %s folder with pattern %s: %s', folderPath, extensions{k}, ME_dir.message));
                set(hPopup, 'String', {['Read Error']}, 'Enable', 'off', 'Value', 1);
                % Optionally return here if one extension fails, or continue
                % return;
            end
        end
    end

    fileList = unique(fileList); % Remove duplicates

    % --- Update Popup Menu ---
    placeholder = ['(Select ', typeName, ')'];
    if isempty(fileList)
        displayList = {['No ', typeName, ' files found']};
        isEnabled = 'off';
        logStatusHelper(handles_caller, ['Warning: No ', typeName, ' files found in ', folderPath]);
    else
        % Sort list alphabetically (case-insensitive)
        [~, sortIdx] = sort(lower(fileList));
        fileListSorted = fileList(sortIdx);
        % Prepend placeholder
        displayList = [ {placeholder}; fileListSorted(:) ];
        isEnabled = 'on';
        logStatusHelper(handles_caller, sprintf('Found %d %s files.', length(fileList), typeName));
    end

    % Set final popup properties
    set(hPopup, 'String', displayList);
    set(hPopup, 'Value', 1); % Select placeholder (index 1) initially
    set(hPopup, 'Enable', isEnabled);

end % End function populateDropdown


% =========================================================
% --- Log Status Helper (CORRECTED) ---
% =========================================================
function logStatusHelper(handles_caller, message)
% Logs message using the main GUI's log box via the 'figure' field in handles_caller.
% INPUT: handles_caller - Handles struct from the calling GUI (MUST contain 'figure')

    % --- Get Main GUI Handles ---
    mainHandles = []; % Initialize
    hMainFigure = []; % Initialize
    % --- CORRECTED CHECK: Look for 'figure' field ---
    if isstruct(handles_caller) && isfield(handles_caller, 'figure') && ishandle(handles_caller.figure)
        hMainFigure = handles_caller.figure; % Get main figure handle from 'figure' field
        try
            mainHandles = guidata(hMainFigure); % Get main GUI handles struct
        catch ME_guidata
             fprintf(2,'logStatusHelper (within populateSingleDropdown): Failed to get guidata from main figure: %s\n', ME_guidata.message);
             % Continue without mainHandles, will log to console
        end
    else
         % --- UPDATED WARNING MESSAGE ---
         fprintf(2,'logStatusHelper (within populateSingleDropdown): Invalid handles_caller struct or missing/invalid ''figure'' handle provided.\n');
         % Continue without mainHandles, will log to console
    end

    % --- Log Message ---
    logBoxHandle = [];
    if isstruct(mainHandles) && isfield(mainHandles, 'logEditBox')
        if ~isempty(mainHandles.logEditBox) && ishandle(mainHandles.logEditBox)
            logBoxHandle = mainHandles.logEditBox;
        end
    end

    if ~isempty(logBoxHandle)
        % Log to the main GUI's log box
        try
            hLog = logBoxHandle;
            hLogFig = ancestor(hLog, 'figure'); % Get log figure handle
            timestamp = datestr(now, 'HH:MM:SS');
            currentLog = get(hLog, 'String');
            if isempty(currentLog), currentLog = {}; elseif ~iscell(currentLog), currentLog = {currentLog}; end
            if ischar(message), messageLines = strsplit(message, '\n'); elseif iscell(message), messageLines = message; else messageLines = {'Log Err: Invalid msg type'}; end
            newLogEntries = cell(sum(~cellfun('isempty', strtrim(messageLines))), 1);
            count = 0;
            for i = 1:length(messageLines)
                lineTxt = strtrim(messageLines{i});
                if ~isempty(lineTxt), count=count+1; newLogEntries{count} = sprintf('[%s] %s', timestamp, lineTxt); end
            end
            updatedLog = [currentLog; newLogEntries(1:count)];
            maxLines = 200;
            if length(updatedLog) > maxLines, updatedLog = updatedLog(end-maxLines+1:end); end
            set(hLog, 'String', updatedLog);
            set(hLog, 'Value', length(updatedLog)); % Scroll to bottom
            % Ensure log window is visible
            if ishandle(hLogFig) && strcmp(get(hLogFig,'Visible'), 'off'), set(hLogFig,'Visible','on'); end
            drawnow limitrate; % Update display
        catch ME_log
            fprintf(2, 'Internal logStatusHelper (in populateSingleDropdown) error writing to log box: %s\n', ME_log.message);
            disp(['LOG (Write Error): ', message]); % Fallback to console
        end
    else
        % Fallback to console if main log box handle not found or invalid
        disp(['LOG (No Win/Main Handles): ', message]);
    end
end % End logStatusHelper function