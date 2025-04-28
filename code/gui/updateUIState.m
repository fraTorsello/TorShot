% gui/updateUIState.m
function updateUIState(hFigure)
% Updates the enabled/disabled state of GUI controls based on application state.
% INPUT: hFigure - Handle to the main simulator figure.

    % --- Input Validation ---
    if nargin < 1 || ~ishandle(hFigure)
        warning('updateUIState: Invalid main figure handle provided.');
        return;
    end
    try
        handles = guidata(hFigure);
        if isempty(handles)
            warning('updateUIState: Could not retrieve guidata from figure.');
            return;
        end
    catch ME_guidata
        warning('updateUIState: Error getting guidata: %s', ME_guidata.message);
        return;
    end

    % --- Determine Current State ---
    componentsLoaded = isfield(handles, 'isComponentsLoaded') && handles.isComponentsLoaded;
    % Check if barrelData exists and is a non-empty struct
    barrelLoaded = isfield(handles, 'barrelData') && isstruct(handles.barrelData) && ~isempty(fieldnames(handles.barrelData));
    paramsCalculated = isfield(handles, 'parametri') && ~isempty(handles.parametri);

    fprintf('DEBUG (updateUIState): componentsLoaded=%d, barrelLoaded=%d, paramsCalculated=%d\n', componentsLoaded, barrelLoaded, paramsCalculated); % Debug print

    % --- Enable/Disable Controls ---

    % Component Selection (Always enabled to allow re-selection)
    trySetEnable(handles, 'pushbuttonLoadComponents', 'on');
    trySetEnable(handles, 'popupmenuCartridge', 'on');
    trySetEnable(handles, 'popupmenuBullet', 'on');
    trySetEnable(handles, 'popupmenuPowder', 'on');
    trySetEnable(handles, 'pushbuttonSelectBarrel', 'on'); % Barrel selection always possible

    % Editable Fields (Enabled only if base components are loaded)
    editableFields = { 'editBulletName', 'editBulletMass_gr', 'editBulletLength_in', 'editBulletDiameter_in', ...
                       'editCartridgeName', 'editCaseLength_in', 'editMaxCaseCapacity_grH2O', 'editBoreDiameter_in', ...
                       'editSeatingDepth_in', ...
                       'editPropellantCharge_gr', 'editBarrelLength_in', 'editTwistRate_inPerTurn', ...
                       'editInitialTemperature_K', 'editAmbientPressure_Pa', 'editShotStartPressure_MPa' };
    stateLoad = 'off';
    bgColor = [0.94 0.94 0.94]; % Disabled background color
    if componentsLoaded
        stateLoad = 'on';
        bgColor = [1 1 1]; % Enabled background color
    end
    for i = 1:length(editableFields)
        if isfield(handles, editableFields{i}) && ishandle(handles.(editableFields{i}))
            set(handles.(editableFields{i}), 'Enable', stateLoad);
            set(handles.(editableFields{i}), 'BackgroundColor', bgColor);
        end
    end

    % Apply & Calc Button (Enabled only if base components AND barrel are loaded)
    stateApply = 'off';
    if componentsLoaded && barrelLoaded
        stateApply = 'on';
    end
    trySetEnable(handles, 'pushbuttonApplyAndCalc', stateApply);
    fprintf('DEBUG (updateUIState): Apply Button State: %s\n', stateApply); % Debug print

    % Run Simulation Button (Enabled only if parameters have been calculated)
    stateRun = 'off';
    if paramsCalculated
        stateRun = 'on';
    end
    trySetEnable(handles, 'pushbuttonRunSimulation', stateRun);

    % Powder Parameters Button/Window
    powderButtonExists = isfield(handles, 'pushbuttonPowderParams') && ishandle(handles.pushbuttonPowderParams);
    powderWindowExists = isfield(handles, 'powderFigure') && ishandle(handles.powderFigure);
    if powderButtonExists
        trySetEnable(handles, 'pushbuttonPowderParams', stateLoad); % Enable if base components loaded
    end
    if ~componentsLoaded && powderWindowExists % Hide powder window if components not loaded
        try set(handles.powderFigure, 'Visible', 'off'); catch; end
    end

    drawnow; % Ensure UI updates are displayed

end % End function updateUIState

% --- Helper function to safely set Enable property ---
function trySetEnable(handles, controlTag, state)
    if isfield(handles, controlTag) && ishandle(handles.(controlTag))
        try
            set(handles.(controlTag), 'Enable', state);
        catch ME_set
             fprintf(2, 'Warning: Could not set Enable="%s" for control "%s": %s\n', state, controlTag, ME_set.message);
        end
    % else
    %     fprintf('Debug: Control tag "%s" not found in handles for trySetEnable.\n', controlTag);
    end
end