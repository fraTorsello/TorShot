% gui/logStatus.m
function logStatus(message)
% Scrive un messaggio nella finestra di log separata (se esiste)
% o nel Command Window come fallback.
% INPUT: message - Stringa o cell array di stringhe da loggare.

    try
        % Cerca l'handle della casella di testo del log tramite il suo Tag
        logBoxHandle = findobj('Tag', 'logEditBox');

        if ~isempty(logBoxHandle) && ishandle(logBoxHandle(1)) % Usa il primo trovato se ce n'è più di uno
            hLog = logBoxHandle(1);
            hLogFig = ancestor(hLog, 'figure'); % Ottieni la figura del log

            if ~ishandle(hLogFig)
                disp(['LOG (Log Fig Invalid): ', message]);
                return;
            end

            % Procedi con il logging nella finestra separata
            timestamp = datestr(now, 'HH:MM:SS');
            currentLog = get(hLog, 'String');
            if isempty(currentLog), currentLog = {}; elseif ~iscell(currentLog), currentLog = {currentLog}; end

            if ischar(message), messageLines = strsplit(message, '\n');
            elseif iscell(message), messageLines = message;
            else messageLines = {'Log Err: Invalid msg type'}; end

            newLogEntries = {};
            for i = 1:length(messageLines)
                lineTxt = strtrim(messageLines{i});
                if ~isempty(lineTxt), newLogEntries{end+1,1} = sprintf('[%s] %s', timestamp, lineTxt); end
            end

            updatedLog = [currentLog; newLogEntries];
            maxLines = 200;
            if length(updatedLog) > maxLines, updatedLog = updatedLog(end-maxLines+1:end); end

            set(hLog, 'String', updatedLog);
            set(hLog, 'Value', length(updatedLog)); % Scrolla in fondo

            % Assicura visibilità finestra log
            if strcmp(get(hLogFig,'Visible'), 'off')
                set(hLogFig,'Visible','on');
            end
            drawnow limitrate;
        else
            % Fallback al Command Window se la finestra non è trovata
            disp(['LOG (No Win): ', message]);
        end
    catch ME_log
        fprintf(2, 'Internal logStatus error: %s\n', ME_log.message);
        fprintf(2, 'Original Log Message: %s\n', message);
        disp(['LOG (Error): ', message]); % Mostra comunque il messaggio
    end
end