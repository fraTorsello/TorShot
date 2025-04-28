% gui/clearPlots.m
function clearPlots(hFigure)
% Pulisce tutti gli assi dei grafici dei risultati nella GUI principale.
% INPUT: hFigure - Handle della figura principale della GUI.

    handles = guidata(hFigure);
    % logStatus(hFigure, 'Clearing plots...'); % Errato, cambia in:
    logStatus('Clearing plots...'); % Corretto

    % ... (resto del codice per definire axesHandles, axesTags, etc.) ...
    axesHandles = {handles.axesPressure, handles.axesPosition, handles.axesVelocity, handles.axesOmega, handles.axesMass, handles.axesEnergy};
    axesTags = {'axesPressure', 'axesPosition', 'axesVelocity', 'axesOmega', 'axesMass', 'axesEnergy'};
    axesTitles = {'Pressure', 'Position', 'Velocity', 'Angular Velocity', 'Masses', 'Energy/Work'};
    axesXLabels = {'Time [ms]', 'Time [ms]', 'Time [ms]', 'Time [ms]', 'Time [ms]', 'Time [ms]'};
    axesYLabels = {'Pressure [MPa]', 'Position [mm]', 'Velocity [m/s]', 'Velocity [RPM]', 'Mass [g]', 'Energy [kJ]'};

    for i = 1:length(axesHandles)
        ax = axesHandles{i};
        axTag = axesTags{i};
        if isfield(handles, axTag) && ishandle(ax)
            try
                cla(ax, 'reset'); title(ax, axesTitles{i}); xlabel(ax, axesXLabels{i}); ylabel(ax, axesYLabels{i}); grid(ax, 'on'); legend(ax, 'off');
            catch ME_cla
                % logStatus(hFigure, sprintf('Warn: Could not clear axes "%s": %s', axTag, ME_cla.message)); % Errato, cambia in:
                logStatus(sprintf('Warn: Could not clear axes "%s": %s', axTag, ME_cla.message)); % Corretto
            end
        else
             % logStatus(hFigure, sprintf('Warn: Axes handle "%s" not found/invalid during clear.', axTag)); % Errato, cambia in:
             logStatus(sprintf('Warn: Axes handle "%s" not found/invalid during clear.', axTag)); % Corretto
        end
    end
    % logStatus(hFigure, 'Plots cleared.'); % Errato, cambia in:
    logStatus('Plots cleared.'); % Corretto
end