function figureCloseRequest(hObject, eventdata)
% Callback executes when the user attempts to close the figure window.

hFigure = hObject; % For close request, hObject is the figure handle
try
    % Attempt to log closure using the external logStatus function
    handles = guidata(hFigure); % Get handles if available
    if ~isempty(handles) && isfield(handles, 'editStatusLog') && ishandle(handles.editStatusLog)
         logStatus(hFigure, 'Chiusura GUI.'); % CALLS EXTERNAL logStatus
    else
        disp('Chiusura GUI (log non disponibile).');
    end
catch ME
    disp('Chiusura GUI (errore durante il log).');
    fprintf(2,'Log Error: %s\n', ME.message);
end
delete(hFigure); % Close the figure window
end