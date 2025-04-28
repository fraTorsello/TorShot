% gui/populateComponentDropdowns.m
function populateComponentDropdowns(hFigure)
% Scansiona le cartelle di configurazione e popola i menu dropdown
% per Cartridge, Bullet, e Powder.
% CORREZIONE: Aggiunge la riga placeholder "(Select...)" all'inizio.
% CORREZIONE: Calcola baseDir salendo dalla cartella 'gui'.
% INPUT: hFigure - Handle della figura principale della GUI.

    handles = guidata(hFigure);
    localLogStatus(handles, 'Scanning for components...');

    % --- Define Folders ---
    try
        % Get base directory robustly
        if isdeployed
            % Logic for deployed application (might need adjustment based on compiler/OS)
            [status, sysPath] = system('path');
             if status == 0 && ~isempty(sysPath)
                 pathParts = strsplit(sysPath, pathsep);
                 baseDir = fileparts(pathParts{1});
             else
                 baseDir = pwd; % Fallback
             end
             if isempty(baseDir) || strcmp(baseDir,'.'), baseDir = pwd; end
             fprintf('DEBUG: Detected deployed application. Base directory estimated as: %s\n', baseDir);
        else
             % Logic for running from MATLAB editor/command window
             try
                 scriptPath = mfilename('fullpath');
                 if isempty(scriptPath)
                     baseDir = pwd;
                     fprintf('DEBUG: mfilename(''fullpath'') returned empty. Using pwd: %s\n', baseDir);
                 else
                     scriptDir = fileparts(scriptPath); % Gets the 'gui' directory
                     % *** FIX: Go UP one level from the script's directory ***
                     baseDir = fileparts(scriptDir);
                     fprintf('DEBUG: mfilename(''fullpath'') used. Script dir: %s. Base dir (parent): %s\n', scriptDir, baseDir);
                 end
             catch ME_mfilename
                 baseDir = pwd; % Fallback
                 fprintf('DEBUG: Error getting script path (%s). Using pwd: %s\n', ME_mfilename.message, baseDir);
             end
             if isempty(baseDir) || strcmp(baseDir,'.'), baseDir = pwd; end
        end
        fprintf('DEBUG: Final baseDir used for config paths: %s\n', baseDir);

        % --- Define full paths based on the CORRECTED baseDir ---
        cartridgeFolder = fullfile(baseDir, 'config', 'cartridges');
        bulletFolder    = fullfile(baseDir, 'config', 'bullets');
        powderFolder    = fullfile(baseDir, 'config', 'powders');

        fprintf('DEBUG: Attempting to populate Cartridge from: %s\n', cartridgeFolder);
        fprintf('DEBUG: Attempting to populate Bullet from: %s\n', bulletFolder);
        fprintf('DEBUG: Attempting to populate Powder from: %s\n', powderFolder);

    catch ME_pathfind
        localLogStatus(handles, ['ERROR finding config folders: ', ME_pathfind.message]);
        fprintf(2, 'CRITICAL ERROR determining base directory or config paths:\n%s\n', ME_pathfind.getReport);
        return; % Cannot proceed without folders
    end

    % --- Populate Dropdowns ---
    populateDropdown(handles.popupmenuCartridge, cartridgeFolder, {'*.m'}, 'Cartridge', handles);
    populateDropdown(handles.popupmenuBullet,    bulletFolder,    {'*.m'}, 'Bullet', handles);
    populateDropdown(handles.popupmenuPowder,    powderFolder,    {'*.mat'}, 'Powder', handles);

    localLogStatus(handles, 'Component scan complete.');

end % Fine populateComponentDropdowns


% =========================================================
% --- Funzione interna (o esterna in utils/) per popolare un singolo dropdown ---
% =========================================================
function populateDropdown(hPopup, folder, extensions, typeName, handles)
% Populates a single popup menu by scanning a folder.
% (Codice interno di questa funzione rimane invariato rispetto alla versione precedente)

    set(hPopup, 'String', {'Scanning...'}, 'Value', 1, 'Enable', 'off'); drawnow;
    fileList = {};
    foundFolder = false;
    folderPath = folder; % Usa il percorso già calcolato

    fprintf('DEBUG (populateDropdown %s): Checking provided path: %s\n', typeName, folderPath);

    % --- Check Folder Existence ---
    if exist(folderPath, 'dir')
        foundFolder = true;
        fprintf('DEBUG (populateDropdown %s): Folder exists: %s\n', typeName, folderPath);
    else
        localLogStatus(handles, sprintf('ERROR: %s folder NOT FOUND at calculated path: %s', typeName, folderPath));
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
                 fprintf('DEBUG (populateDropdown %s): Found %d items matching "%s" in %s\n', typeName, length(files), currentExtension, folderPath);
                for i = 1:length(files)
                    if files(i).isdir, continue; end
                    [~, name, ~] = fileparts(files(i).name);
                    if ~isempty(name) && ~startsWith(name, '.') && ~endsWith(name, '~')
                         fileList{end+1} = name; %#ok<AGROW>
                    end
                end
            catch ME_dir
                localLogStatus(handles, sprintf('Warning: Could not read %s folder with pattern %s: %s', folderPath, extensions{k}, ME_dir.message));
                set(hPopup, 'String', {['Read Error']}, 'Enable', 'off', 'Value', 1);
            end
        end
    end

    fileList = unique(fileList);

    % --- Update Popup Menu ---
    placeholder = ['(Select ', typeName, ')'];
    if isempty(fileList)
        displayList = {['No ', typeName, ' files found']};
        isEnabled = 'off';
        localLogStatus(handles, ['Warning: No ', typeName, ' files found in ', folderPath]);
    else
        [~, sortIdx] = sort(lower(fileList));
        fileListSorted = fileList(sortIdx);
        displayList = [ {placeholder}; fileListSorted(:) ];
        isEnabled = 'on';
        localLogStatus(handles, sprintf('Found %d %s files.', length(fileList), typeName));
    end

    set(hPopup, 'String', displayList);
    set(hPopup, 'Value', 1);
    set(hPopup, 'Enable', isEnabled);

end

% --- Local Log Status Helper ---
function localLogStatus(handles, message)
    if isfield(handles, 'logEditBox') && ~isempty(handles.logEditBox) && ishandle(handles.logEditBox)
        try
            hLog = handles.logEditBox;
            hLogFig = ancestor(hLog, 'figure');
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
            set(hLog, 'Value', length(updatedLog));
            if ishandle(hLogFig) && strcmp(get(hLogFig,'Visible'), 'off'), set(hLogFig,'Visible','on'); end
            drawnow limitrate;
        catch ME_log
            fprintf(2, 'Internal localLogStatus error: %s\n', ME_log.message);
            disp(['LOG (Error): ', message]);
        end
    else
        disp(['LOG (No Win): ', message]);
    end
end
