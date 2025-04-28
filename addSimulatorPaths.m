% gui/addSimulatorPaths.m
function addSimulatorPaths()
% Aggiunge le cartelle necessarie del simulatore al path di MATLAB.

    try
        % Tenta di determinare la directory dello script corrente o CWD
        if isdeployed % Codice compilato
            [~, sysPath] = system('path'); % Ottiene il path dell'eseguibile (potrebbe richiedere aggiustamenti)
             baseDir = fileparts(sysPath);
             if isempty(baseDir) || strcmp(baseDir,'.'), baseDir = pwd; end
        else % Codice interpretato
             try
                baseDir = fileparts(mfilename('fullpath')); % Directory dello script .m
             catch
                 baseDir = pwd; % Fallback
             end
             if isempty(baseDir) || strcmp(baseDir,'.'), baseDir = pwd; end
        end
        fprintf('Base application directory identified as: %s\n', baseDir);

        pathsToAdd = { ...
            fullfile(baseDir, 'config'), ...
            fullfile(baseDir, 'config', 'cartridges'), ...
            fullfile(baseDir, 'config', 'bullets'), ...
            fullfile(baseDir, 'config', 'powders'), ...
            fullfile(baseDir, 'core'), ...
            fullfile(baseDir, 'utils'), ...
            fullfile(baseDir, 'postprocessing'), ...
            fullfile(baseDir, 'output'), ...
            fullfile(baseDir, 'output', 'Matlab_Simulation_Data'), ...
            fullfile(baseDir, 'gui') ...
        };
        pathsAddedCount = 0;
        currentPath = path; % Get current path string
        pathsToAddFiltered = {}; % Store paths that are not already added

        for i = 1:length(pathsToAdd)
            normalizedPath = pathsToAdd{i};
            % Controlla se il path esiste già nel path di MATLAB
            pathExists = contains(path, [normalizedPath pathsep], 'IgnoreCase', ispc) || ...
                         strcmp(path, normalizedPath) || ... % Caso esatto
                         endsWith(path, [pathsep normalizedPath], 'IgnoreCase', ispc); % Caso finisca con esso

            if ~pathExists
                 pathsToAddFiltered{end+1} = normalizedPath; %#ok<AGROW>
            end
        end

        for i = 1:length(pathsToAddFiltered)
             p = pathsToAddFiltered{i};
             if exist(p, 'dir')
                 addpath(p);
                 pathsAddedCount = pathsAddedCount + 1;
                 fprintf('Added path: %s\n', p);
             else
                 warning('Simulator path not found, cannot add: %s', p);
             end
        end

        if pathsAddedCount > 0
            rehash path; % Refresh path cache only if paths were added
            fprintf('Added %d new simulator paths.\n', pathsAddedCount);
        else
            disp('All required simulator paths already seem to be in the MATLAB path.');
        end
    catch ME
        warning('Error adding simulator paths: %s', ME.message);
         fprintf(2, 'Stack trace for path error:\n'); disp(ME.getReport);
    end
end