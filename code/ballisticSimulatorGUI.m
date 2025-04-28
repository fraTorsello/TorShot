function ballisticSimulatorGUI()
% Internal Ballistics Simulator GUI (Responsive & Dynamic Layout)
% --- REFACTORED for Component-Based Selection & Separate Log Window ---
% --- VERSIONE CON FUNZIONI HELPER ESTERNALIZZATE ---
% --- VERSIONE CON DISPLAY/INPUT IN GRANI/POLLICI E CAMPI EDITABILI ---
% --- Logica OAL Calcolata da Seating Depth ---
% --- FIX: Corrected layout, button logic, added local log fallback ---
% --- MODIFIED: Added call to update barrel stress window ---

%% --- Add necessary paths ---
disp('Adding paths for GUI and Simulator...');
if exist('addSimulatorPaths', 'file') == 2
    addSimulatorPaths();
else
    warning('addSimulatorPaths function not found.');
    % Fallback path addition (adjust as necessary)
    try
        appDir = fileparts(mfilename('fullpath')); if isempty(appDir), appDir=pwd; end
        folders = {'gui', 'utils', 'core', 'postprocessing', 'config', fullfile('config','cartridges'), fullfile('config','bullets'), fullfile('config','powders')};
        for i = 1:length(folders)
            folderPath = fullfile(appDir, folders{i});
            if exist(folderPath,'dir'), addpath(folderPath); fprintf('Manually added path: %s\n', folderPath); end
        end
        rehash path;
    catch ME_path
        warning('Failed to add paths manually: %s', ME_path.message);
    end
end

%% --- Create Main Window ---
handles = struct();
mainFigPosition = [0.05 0.05 0.65 0.85]; % Normalized units
handles.figure = figure('Name', 'Internal Ballistics Simulator (Components)', ...
                        'Units', 'normalized', 'Position', mainFigPosition, ...
                        'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
                        'Resize', 'on', 'ResizeFcn', @onMainFigureResize, ...
                        'CloseRequestFcn', @mainFigureCloseRequest, ...
                        'Visible', 'off', 'Tag', 'mainSimulatorFigure');

% --- Initialize Other GUIs ---
% Powder Parameters GUI
if exist('powderParametersGUI', 'file') == 2
    handles.powderFigure = powderParametersGUI(handles.figure);
else
    handles.powderFigure = [];
    warning('powderParametersGUI.m not found.');
end
% Log Window GUI
if exist('createLogWindow', 'file') == 2
    handles.logEditBox = createLogWindow(mainFigPosition);
    if ishandle(handles.logEditBox)
        handles.logFigure = ancestor(handles.logEditBox, 'figure');
    else
        handles.logFigure = []; % Handle case where log window creation failed
    end
else
    handles.logEditBox = [];
    handles.logFigure = [];
    warning('createLogWindow.m not found.');
end

guidata(handles.figure, handles); % Save initial handles

%% --- Initialize Data ---
handles.componentData = [];
handles.barrelData = []; % Initialize barrel data field
handles.parametri = [];
handles.isComponentsLoaded = false;
handles.isSimRun = false;
handles.stressResults = []; % Initialize stress results field

%% --- Define Constants (Include Inverse) ---
handles.CONV = struct(...
    'KG_PER_GRAIN', 1/15432.35835, ...
    'METERS_PER_INCH', 0.0254, ...
    'MM_PER_M', 1000, ...
    'G_PER_KG', 1000, ...
    'M3_PER_GRH2O', (1/15432.35835)/1000.0, ...
    'CM3_PER_M3', 1e6, ...
    'MM2_PER_M2', 1e6, ...
    'MPA_PER_PA', 1e-6, ...
    'PA_PER_MPA', 1e6, ...
    'GRAINS_PER_KG', 15432.35835, ...
    'INCHES_PER_METER', 1/0.0254, ...
    'GRH2O_PER_M3', 1000.0*15432.35835, ...
    'IN2_PER_M2', (1/0.0254)^2 ...
);

%% --- Create Main Panels ---
handles.leftPanel = uipanel('Parent', handles.figure, 'Title', 'Component Selection & Parameters', ...
                           'Units', 'normalized', 'Position', [0.02 0.02 0.45 0.96], 'Tag', 'leftPanel');
handles.rightPanel = uipanel('Parent', handles.figure, 'Title', 'Simulation Results', ...
                            'FontWeight', 'bold', 'Units', 'normalized', ...
                            'Position', [0.49 0.02 0.50 0.96], 'Tag', 'rightPanel');

%% --- Create Controls in Left Panel ---
panelMargin = 0.03;
currentY = 1.0 - panelMargin; % Start from top

% --- Selection Panel ---
selHeight = 0.18; % Relative height
selYpos = currentY - selHeight;
handles.panelSelection = uipanel('Parent', handles.leftPanel,'Title', 'Select Components', ...
                                 'Units', 'normalized', 'Position', [panelMargin, selYpos, 1-2*panelMargin, selHeight], ...
                                 'Tag', 'panelSelection');
currentY = selYpos - 0.01; % Update Y position below panel

% Controls inside Selection Panel
selMarginX = 0.04; selLabelW = 0.25; selDropW = 0.68; selH = 0.25; % Normalized within panel
selY = 1.0 - selH - 0.05; % Start Y inside panel
uicontrol('Parent', handles.panelSelection, 'Style', 'text', 'String', 'Cartridge:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [selMarginX, selY, selLabelW, selH]);
handles.popupmenuCartridge = uicontrol('Parent', handles.panelSelection, 'Style', 'popupmenu', 'String', {'(Scan...)'}, 'Units', 'normalized', 'Position', [selMarginX+selLabelW, selY, selDropW, selH], 'Tag', 'popupmenuCartridge');
selY = selY - selH - 0.03;
uicontrol('Parent', handles.panelSelection, 'Style', 'text', 'String', 'Bullet:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [selMarginX, selY, selLabelW, selH]);
handles.popupmenuBullet = uicontrol('Parent', handles.panelSelection, 'Style', 'popupmenu', 'String', {'(Scan...)'}, 'Units', 'normalized', 'Position', [selMarginX+selLabelW, selY, selDropW, selH], 'Tag', 'popupmenuBullet');
selY = selY - selH - 0.03;
uicontrol('Parent', handles.panelSelection, 'Style', 'text', 'String', 'Powder:', 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [selMarginX, selY, selLabelW, selH]);
handles.popupmenuPowder = uicontrol('Parent', handles.panelSelection, 'Style', 'popupmenu', 'String', {'(Scan...)'}, 'Units', 'normalized', 'Position', [selMarginX+selLabelW, selY, selDropW, selH], 'Tag', 'popupmenuPowder');

% --- Load Components Button ---
loadButtonHeight = 0.04; % Relative height
loadButtonY = currentY - loadButtonHeight;
handles.pushbuttonLoadComponents = uicontrol('Parent', handles.leftPanel, 'Style', 'pushbutton', ...
                                             'String', 'Load Selected Components', 'Units', 'normalized', ...
                                             'Position', [panelMargin, loadButtonY, 1-2*panelMargin, loadButtonHeight], ...
                                             'Callback', @loadComponentsCallback, 'Tag', 'pushbuttonLoadComponents');
currentY = loadButtonY - 0.01; % Space below button

% --- Select Barrel Button ---
selectBarrelButtonHeight = 0.04; % Relative height
selectBarrelButtonY = currentY - selectBarrelButtonHeight;
handles.pushbuttonSelectBarrel = uicontrol('Parent', handles.leftPanel, 'Style', 'pushbutton', ...
                                           'String', 'Select Barrel...', 'Units', 'normalized', ...
                                           'Position', [panelMargin, selectBarrelButtonY, 1-2*panelMargin, selectBarrelButtonHeight], ...
                                           'Callback', @selectBarrelCallback, 'Tag', 'pushbuttonSelectBarrel', ...
                                           'Enable', 'on'); % Initially enabled
currentY = selectBarrelButtonY - 0.01; % Space below button

% --- Powder Parameters Button ---
powderButtonHeight = 0.04; % Relative height
powderButtonY = currentY - powderButtonHeight;
handles.pushbuttonPowderParams = uicontrol('Parent', handles.leftPanel, 'Style', 'pushbutton', ...
                                           'String', 'Display Powder Parameters', 'Units', 'normalized', ...
                                           'Position', [panelMargin, powderButtonY, 1-2*panelMargin, powderButtonHeight], ...
                                           'Callback', @openPowderParamsCallback_Wrapper, ... % Use Wrapper
                                           'Tag', 'pushbuttonPowderParams', 'Enable', 'off'); % Initially disabled
currentY = powderButtonY - 0.02; % Larger gap

% --- Dimensions Panel ---
% Calculate remaining height for Dim and Sim panels dynamically
bottomOfPowderButton = powderButtonY;
bottomOfMainButtons = panelMargin + 0.01 + 0.04 + 0.01 + 0.04; % RunBtn + Gap + ApplyBtn + Gap
availableHeightForPanels = bottomOfPowderButton - bottomOfMainButtons - 0.01; % Available space
dimPanelRelHeight = 0.65; % Proportion for dimensions panel
simPanelRelHeight = 0.35; % Proportion for sim params panel

dimPanelHeight = max(0.1, availableHeightForPanels * dimPanelRelHeight); % Ensure minimum height
dimPanelY = currentY - dimPanelHeight;
handles.panelDimensions = uipanel('Parent', handles.leftPanel, 'Title', 'Dimensions & Loading Parameters', ...
                                  'Units', 'normalized', 'Position', [panelMargin, dimPanelY, 1-2*panelMargin, dimPanelHeight], ...
                                  'Tag', 'panelDimensions');
currentY = dimPanelY - 0.02; % Gap

% Create fields inside Dimensions Panel
dimMarginX = 0.03; dimLabelW = 0.55; dimValueW = 0.38; dimFieldH = 0.055; dimSmallGap = 0.008; dimSectionGap = 0.02;
dimVY = 1.0 - dimFieldH - dimSmallGap; % Start Y inside panel
[handles.labelBulletName, handles.editBulletName] = createEditableField(handles.panelDimensions, 'Bullet Name/ID', 'BulletName', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 'N/A'); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelBulletMass_gr, handles.editBulletMass_gr] = createEditableField(handles.panelDimensions, 'Bullet Weight (gr)', 'BulletMass_gr', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 95.0); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelBulletLength_in, handles.editBulletLength_in] = createEditableField(handles.panelDimensions, 'Bullet Length (in)', 'BulletLength_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 1.031); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelBulletDiam_in, handles.editBulletDiameter_in] = createEditableField(handles.panelDimensions, 'Bullet Diameter (in)', 'BulletDiameter_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 0.243); dimVY = dimVY - dimSectionGap;
[handles.labelCartridgeName, handles.editCartridgeName] = createEditableField(handles.panelDimensions, 'Cartridge Name', 'CartridgeName', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 'N/A'); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelCaseLength_in, handles.editCaseLength_in] = createEditableField(handles.panelDimensions, 'Case Length (in)', 'CaseLength_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 2.044); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelMaxCap_grH2O, handles.editMaxCaseCapacity_grH2O] = createEditableField(handles.panelDimensions, 'Max Case Cap (gr H2O)', 'MaxCaseCapacity_grH2O', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 54.0); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelBoreDiam_in, handles.editBoreDiameter_in] = createEditableField(handles.panelDimensions, 'Bore/Groove Dia. (in)', 'BoreDiameter_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 0.243); dimVY = dimVY - dimSectionGap;
[handles.labelSeatingDepth_in, handles.editSeatingDepth_in] = createEditableField(handles.panelDimensions, 'Seating Depth (in)', 'SeatingDepth_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW, 0.400); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelCartridgeOAL_in, handles.textCartridgeOAL_in] = createReadOnlyField(handles.panelDimensions, 'Cartridge OAL (in)', 'CartridgeOAL_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW); dimVY = dimVY - dimSectionGap;
[handles.labelVolOccBullet_grH2O, handles.textVolumeOccupiedBullet_grH2O] = createReadOnlyField(handles.panelDimensions, 'Vol Occupied Bullet (gr H2O)', 'VolumeOccupiedBullet_grH2O', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelUseableCap_grH2O, handles.textQLUseableCaseCapacity_grH2O] = createReadOnlyField(handles.panelDimensions, 'Useable Case Cap (QL) (gr H2O)', 'QLUseableCaseCapacity_grH2O', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelInitialFreeVolumeV0_grH2O, handles.textInitialFreeVolumeV0_grH2O] = createReadOnlyField(handles.panelDimensions, 'Initial Free Vol V0 (gr H2O)', 'InitialFreeVolumeV0_grH2O', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW); dimVY = dimVY - dimFieldH - dimSmallGap;
[handles.labelBulletTravel_in, handles.textBulletTravel_in] = createReadOnlyField(handles.panelDimensions, 'Bullet Travel (in)', 'BulletTravel_in', dimVY, dimFieldH, dimMarginX, dimLabelW, dimValueW);

% --- Simulation Parameters Panel ---
simPanelHeight = max(0.1, availableHeightForPanels * simPanelRelHeight); % Ensure minimum height
simPanelY = currentY - simPanelHeight;
handles.panelSimParams = uipanel('Parent', handles.leftPanel,'Title', 'Simulation Parameters', ...
                                 'Units', 'normalized', 'Position', [panelMargin, simPanelY, 1-2*panelMargin, simPanelHeight], ...
                                 'Tag', 'panelSimParams');
currentY = simPanelY - 0.02; % Update Y below panel

% Create fields inside Simulation Parameters Panel
spMarginX = 0.04; spLabelW = 0.55; spValueW = 0.38; spNumFields = 6;
spAvailableHeight = 1.0 - 0.06; % Inner margin top/bottom
spTotalGap = (spNumFields - 1) * 0.01;
spH = max(0.05, (spAvailableHeight - spTotalGap) / spNumFields); % Ensure minimum height
spY = 1.0 - 0.03 - spH; % Start Y inside panel
[handles.labelPropCharge_gr, handles.editPropellantCharge_gr] = createEditableField(handles.panelSimParams, 'Propellant Charge (gr)', 'PropellantCharge_gr', spY, spH, spMarginX, spLabelW, spValueW, 40.0); spY = spY - spH - 0.01;
[handles.labelBarrelLength_in, handles.editBarrelLength_in] = createEditableField(handles.panelSimParams, 'Barrel Length (in)', 'BarrelLength_in', spY, spH, spMarginX, spLabelW, spValueW, 26.0); spY = spY - spH - 0.01;
[handles.labelTwistRate_in, handles.editTwistRate_inPerTurn] = createEditableField(handles.panelSimParams, 'Twist Rate (in/turn)', 'TwistRate_inPerTurn', spY, spH, spMarginX, spLabelW, spValueW, 10.0); spY = spY - spH - 0.01;
[handles.labelInitialTemp_K, handles.editInitialTemperature_K] = createEditableField(handles.panelSimParams, 'Initial Temp (K)', 'InitialTemperature_K', spY, spH, spMarginX, spLabelW, spValueW, 293.15); spY = spY - spH - 0.01;
[handles.labelAmbientP_Pa, handles.editAmbientPressure_Pa] = createEditableField(handles.panelSimParams, 'Ambient Pressure (Pa)', 'AmbientPressure_Pa', spY, spH, spMarginX, spLabelW, spValueW, 101325); spY = spY - spH - 0.01;
[handles.labelShotStartPressure_MPa, handles.editShotStartPressure_MPa] = createEditableField(handles.panelSimParams, 'Shot Start Press (MPa)', 'ShotStartPressure_MPa', spY, spH, spMarginX, spLabelW, spValueW, 25.0);

% --- Bottom Buttons ---
buttonHeight = 0.04; % Relative height
buttonY_Run = panelMargin + 0.01;
buttonY_Apply = buttonY_Run + buttonHeight + 0.01;
handles.pushbuttonApplyAndCalc = uicontrol('Parent', handles.leftPanel, 'Style', 'pushbutton', ...
                                           'String', 'Apply & Calculate Volumes', 'Units', 'normalized', ...
                                           'Position', [panelMargin, buttonY_Apply, 1-2*panelMargin, buttonHeight], ...
                                           'Callback', @applyAndCalcCallback_Wrapper, 'Tag', 'pushbuttonApplyAndCalc', ...
                                           'Enable', 'off'); % Initially disabled
handles.pushbuttonRunSimulation = uicontrol('Parent', handles.leftPanel, 'Style', 'pushbutton', ...
                                            'String', 'Run Simulation', 'FontWeight', 'bold', 'Units', 'normalized', ...
                                            'Position', [panelMargin, buttonY_Run, 1-2*panelMargin, buttonHeight], ...
                                            'Callback', @runSimulationCallback, 'Tag', 'pushbuttonRunSimulation', ...
                                            'Enable', 'off'); % Initially disabled

%% --- Create Axes in Right Panel ---
% Use 4x2 layout to accommodate stress plots
axMarginX = 0.05; axMarginY = 0.05; % Adjusted margins for 4 rows
axWidth = (1 - 3*axMarginX) / 2;
axHeight = (1 - 5*axMarginY) / 4; % Height for 4 rows
xPos1 = axMarginX; xPos2 = 2*axMarginX + axWidth;
yPos1 = 1 - axMarginY - axHeight;
yPos2 = yPos1 - axMarginY - axHeight;
yPos3 = yPos2 - axMarginY - axHeight;
yPos4 = yPos3 - axMarginY - axHeight; % Bottom row Y position

handles.axesPressure = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos1, yPos1, axWidth, axHeight], 'Tag', 'axesPressure'); title(handles.axesPressure, 'Pressure'); xlabel(handles.axesPressure, 'Time [ms]'); ylabel(handles.axesPressure, 'Pressure [MPa]'); grid(handles.axesPressure, 'on');
handles.axesVelocity = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos2, yPos1, axWidth, axHeight], 'Tag', 'axesVelocity'); title(handles.axesVelocity, 'Velocity'); xlabel(handles.axesVelocity, 'Time [ms]'); ylabel(handles.axesVelocity, 'Velocity [m/s]'); grid(handles.axesVelocity, 'on'); % Changed units
handles.axesPosition = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos1, yPos2, axWidth, axHeight], 'Tag', 'axesPosition'); title(handles.axesPosition, 'Position'); xlabel(handles.axesPosition, 'Time [ms]'); ylabel(handles.axesPosition, 'Position [cm]'); grid(handles.axesPosition, 'on'); % Changed units
handles.axesOmega = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos2, yPos2, axWidth, axHeight], 'Tag', 'axesOmega'); title(handles.axesOmega, 'Angular Velocity'); xlabel(handles.axesOmega, 'Time [ms]'); ylabel(handles.axesOmega, 'Velocity [RPM]'); grid(handles.axesOmega, 'on');
handles.axesMass = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos1, yPos3, axWidth, axHeight], 'Tag', 'axesMass'); title(handles.axesMass, 'Masses'); xlabel(handles.axesMass, 'Time [ms]'); ylabel(handles.axesMass, 'Mass [g]'); grid(handles.axesMass, 'on');
handles.axesEnergy = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos2, yPos3, axWidth, axHeight], 'Tag', 'axesEnergy'); title(handles.axesEnergy, 'Energy/Work'); xlabel(handles.axesEnergy, 'Time [ms]'); ylabel(handles.axesEnergy, 'Energy [kJ]'); grid(handles.axesEnergy, 'on'); % Changed units
handles.axesStress = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos1, yPos4, axWidth, axHeight], 'Tag', 'axesStress'); title(handles.axesStress, 'Barrel Stress (Inner)'); xlabel(handles.axesStress, 'Time [ms]'); ylabel(handles.axesStress, 'Stress [MPa]'); grid(handles.axesStress, 'on');
handles.axesSafetyFactor = axes('Parent', handles.rightPanel, 'Units','normalized', 'Position', [xPos2, yPos4, axWidth, axHeight], 'Tag', 'axesSafetyFactor'); title(handles.axesSafetyFactor, 'Safety Factor (Inner)'); xlabel(handles.axesSafetyFactor, 'Time [ms]'); ylabel(handles.axesSafetyFactor, 'SF [-]'); grid(handles.axesSafetyFactor, 'on');


%% --- Final Setup ---
guidata(handles.figure, handles); % Save updated handles
try
    if exist('populateComponentDropdowns','file')==2
        populateComponentDropdowns(handles.figure);
    else
        error('Function populateComponentDropdowns not found.');
    end
    updateUIState(handles.figure); % Call initial UI state update
    localLogStatus(handles,'GUI started. Select components and press "Load".'); % Use local log
catch initME
    fprintf(2,'ERR initializing GUI: %s\n', initME.message);
    fprintf(2,'Stack trace:\n'); disp(initME.getReport);
    localLogStatus(handles, ['ERR init GUI: ' initME.message]); % Use local log
    uiwait(msgbox(sprintf('Error initializing GUI:\n%s\nCheck console.', initME.message), 'Init Error', 'error'));
end
set(handles.figure, 'Visible', 'on'); % Make the main figure visible

end % END FUNCTION ballisticSimulatorGUI

%% --- GUI HELPER FUNCTIONS ---

function [hLabel, hValue] = createReadOnlyField(parent, labelStr, tagSuffix, yPos, height, marginX, labelW, valueW)
    hLabel = uicontrol('Parent', parent, 'Style', 'text', 'String', [labelStr ':'], 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [marginX, yPos, labelW, height], 'Tag', ['label' tagSuffix]);
    hValue = uicontrol('Parent', parent, 'Style', 'text', 'String', '---', 'HorizontalAlignment', 'left', 'Units', 'normalized', 'BackgroundColor', [0.9 0.9 0.9], 'Position', [marginX+labelW+0.01, yPos, valueW, height], 'Tag', ['text' tagSuffix]);
end

function [hLabel, hEdit] = createEditableField(parent, labelStr, tagSuffix, yPos, height, marginX, labelW, valueW, defaultValue)
    hLabel = uicontrol('Parent', parent, 'Style', 'text', 'String', [labelStr ':'], 'HorizontalAlignment', 'right', 'Units', 'normalized', 'Position', [marginX, yPos, labelW, height], 'Tag', ['label' tagSuffix]);
    hEdit = uicontrol('Parent', parent, 'Style', 'edit', 'String', num2str(defaultValue), 'HorizontalAlignment', 'left', 'Units', 'normalized', 'BackgroundColor', [1 1 1], 'Position', [marginX+labelW+0.01, yPos, valueW, height], 'Tag', ['edit' tagSuffix], 'Enable', 'off');
end



function mainFigureCloseRequest(hObject, ~)
% Handles closing the main GUI window.
     handles = guidata(hObject);
     % Close associated windows if they exist and are valid handles
     if isfield(handles, 'logFigure') && ~isempty(handles.logFigure) && ishandle(handles.logFigure)
         try delete(handles.logFigure); catch, end
     end
     if isfield(handles, 'powderFigure') && ~isempty(handles.powderFigure) && ishandle(handles.powderFigure)
         try delete(handles.powderFigure); catch, end
     end
     % Also close the barrel selection window if it exists
     if isfield(handles, 'barrelSelectionFigure') && ~isempty(handles.barrelSelectionFigure) && ishandle(handles.barrelSelectionFigure)
         try delete(handles.barrelSelectionFigure); catch, end
     end
     delete(hObject); % Delete the main figure
     disp('Simulator GUI closed.');
 end

function onMainFigureResize(hObject, ~)
% Callback for main figure resize event.
    try
        updateParameterLayout(hObject); % Update left panel layout
        % Add code here if right panel layout needs dynamic adjustment too
        % Example: Adjust right panel axes based on new figure size
        handles = guidata(hObject);
        if isfield(handles,'rightPanel') && ishandle(handles.rightPanel)
            % Redefine axes positions based on current rightPanel size
            axMarginX = 0.05; axMarginY = 0.05;
            axWidth = (1 - 3*axMarginX) / 2;
            axHeight = (1 - 5*axMarginY) / 4; % 4 rows
            xPos1 = axMarginX; xPos2 = 2*axMarginX + axWidth;
            yPos1 = 1 - axMarginY - axHeight;
            yPos2 = yPos1 - axMarginY - axHeight;
            yPos3 = yPos2 - axMarginY - axHeight;
            yPos4 = yPos3 - axMarginY - axHeight;

            set(handles.axesPressure, 'Position', [xPos1, yPos1, axWidth, axHeight]);
            set(handles.axesVelocity, 'Position', [xPos2, yPos1, axWidth, axHeight]);
            set(handles.axesPosition, 'Position', [xPos1, yPos2, axWidth, axHeight]);
            set(handles.axesOmega, 'Position', [xPos2, yPos2, axWidth, axHeight]);
            set(handles.axesMass, 'Position', [xPos1, yPos3, axWidth, axHeight]);
            set(handles.axesEnergy, 'Position', [xPos2, yPos3, axWidth, axHeight]);
            if isfield(handles,'axesStress') && ishandle(handles.axesStress)
                 set(handles.axesStress, 'Position', [xPos1, yPos4, axWidth, axHeight]);
            end
            if isfield(handles,'axesSafetyFactor') && ishandle(handles.axesSafetyFactor)
                 set(handles.axesSafetyFactor, 'Position', [xPos2, yPos4, axWidth, axHeight]);
            end
        end
    catch ME_resize
        fprintf('Error during layout update on resize: %s\n', ME_resize.message);
    end
end

function updateParameterLayout(hFigure)
% Dynamically redraws the layout of the left panel controls.
    handles = guidata(hFigure);
    if isempty(handles) || ~isfield(handles,'leftPanel') || ~ishandle(handles.leftPanel), return; end

    panelMargin = 0.03;
    buttonHeightRel = 0.04; % Relative height for buttons
    buttonGapRel = 0.01;    % Relative gap between buttons

    % Calculate Y positions for bottom buttons (Apply & Run)
    buttonY_Run = panelMargin + buttonGapRel;
    buttonY_Apply = buttonY_Run + buttonHeightRel + buttonGapRel;
    topOfBottomButtonArea = buttonY_Apply + buttonHeightRel + buttonGapRel; % Y coord where panels above should end

    % Set positions for bottom buttons
    set(handles.pushbuttonRunSimulation, 'Position', [panelMargin, buttonY_Run, 1-2*panelMargin, buttonHeightRel]);
    set(handles.pushbuttonApplyAndCalc, 'Position', [panelMargin, buttonY_Apply, 1-2*panelMargin, buttonHeightRel]);

    % Calculate total available height for panels and buttons above Apply/Run
    topBoundary = 1.0 - panelMargin;
    bottomBoundaryForPanels = topOfBottomButtonArea;
    totalAvailableHeight = topBoundary - bottomBoundaryForPanels;

    % Define controls to manage visibility and layout
    % Order matters for sequential positioning!
    controlsToLayout = {'panelSelection', 'pushbuttonLoadComponents', 'pushbuttonSelectBarrel', 'pushbuttonPowderParams', 'panelDimensions', 'panelSimParams'};
    % Define relative heights (these should sum close to 1, representing proportions)
    relativeHeights = [0.16, 0.04, 0.04, 0.04, 0.42, 0.30]; % Adjusted proportions

    % Check for minimum space
    minRequiredHeight = 0.2; % Minimum total height to display everything
    isVisible = totalAvailableHeight >= minRequiredHeight;

    % Toggle visibility based on available space
    for i = 1:length(controlsToLayout)
        if isfield(handles, controlsToLayout{i}) && ishandle(handles.(controlsToLayout{i}))
             set(handles.(controlsToLayout{i}), 'Visible', iff(isVisible, 'on', 'off'));
        end
    end

    if ~isVisible
        warning('Insufficient vertical space in the left panel to display all controls.');
        return; % Don't proceed with layout if too small
    end

    % Calculate actual heights and gaps
    numElements = length(controlsToLayout);
    numGaps = numElements - 1;
    panelGapRel = 0.01; % Relative gap between panels/buttons

    % Normalize relative heights
    totalRelHeightFactor = sum(relativeHeights);
    actualHeights = (totalAvailableHeight - numGaps * panelGapRel) * (relativeHeights / totalRelHeightFactor);

    % Ensure heights are not negative
    actualHeights(actualHeights < 0.01) = 0.01; % Set a minimum height

    % Adjust heights proportionally if they exceed available space after setting minimums
    currentTotalHeight = sum(actualHeights) + numGaps * panelGapRel;
    if currentTotalHeight > totalAvailableHeight
        scaleDown = totalAvailableHeight / currentTotalHeight;
        actualHeights = actualHeights * scaleDown;
        panelGapRel = panelGapRel * scaleDown;
    end

    % --- POSITION ELEMENTS SEQUENTIALLY from top ---
    currentTopY = topBoundary;
    for i = 1:length(controlsToLayout)
        controlTag = controlsToLayout{i};
        controlHeight = actualHeights(i);
        currentBottomY = currentTopY - controlHeight;

        if isfield(handles, controlTag) && ishandle(handles.(controlTag))
            set(handles.(controlTag), 'Position', [panelMargin, currentBottomY, 1-2*panelMargin, controlHeight]);
            % Update inner layout for panels
            if startsWith(controlTag, 'panel') % && exist('updatePanelInnerLayout','file')==2 % Removed check for simplicity
                if strcmp(controlTag, 'panelDimensions')
                    updatePanelInnerLayout(handles.(controlTag), {'BulletName', 'BulletMass_gr', 'BulletLength_in', 'BulletDiameter_in', 'GAP', 'CartridgeName', 'CaseLength_in', 'MaxCaseCapacity_grH2O', 'BoreDiameter_in', 'GAP', 'SeatingDepth_in', 'CartridgeOAL_in', 'GAP', 'VolumeOccupiedBullet_grH2O', 'QLUseableCaseCapacity_grH2O', 'InitialFreeVolumeV0_grH2O', 'BulletTravel_in'});
                elseif strcmp(controlTag, 'panelSimParams')
                     updatePanelInnerLayout(handles.(controlTag), {'PropellantCharge_gr', 'BarrelLength_in', 'TwistRate_inPerTurn', 'InitialTemperature_K', 'AmbientPressure_Pa', 'ShotStartPressure_MPa'});
                % Add inner layout update for panelSelection if needed
                elseif strcmp(controlTag, 'panelSelection')
                    % Example: Update layout inside selection panel if needed
                    % updateSelectionPanelLayout(handles.panelSelection);
                end
            end
        end
        currentTopY = currentBottomY - panelGapRel; % Move Y for the next element
    end
end

function updatePanelInnerLayout(hPanel, fieldTags)
% Dynamically redraws the layout of elements INSIDE a given panel.
    handles = guidata(hPanel); if isempty(handles), return; end
    try panelPosPixels = getpixelposition(hPanel); panelHeightPixels = panelPosPixels(4);
    catch, panelPosNorm = get(hPanel, 'Position'); panelHeightPixels = panelPosNorm(4) * get(ancestor(hPanel,'figure'), 'Position') * [0;0;0;1]; end % Estimate pixel height

    innerMarginX = 0.03; innerTopMargin = 0.08; innerBottomMargin = 0.05;
    labelW = 0.55; valueW = 0.38; gapX = 0.01;
    minFieldHeightPixels = 17; % Minimum pixel height for a field

    numFields = 0; numGaps = 0;
    for i = 1:length(fieldTags), if ~strcmpi(fieldTags{i}, 'GAP'), numFields = numFields + 1; else numGaps = numGaps + 1; end; end
    if numFields == 0, return; end % No fields to layout

    gapHeightFactor = 0.3; % Gaps are smaller than fields
    totalVerticalUnits = numFields + (numGaps * gapHeightFactor);
    availableHeightNorm = max(0.01, 1.0 - innerTopMargin - innerBottomMargin); % Ensure positive height

    % Calculate ideal heights
    heightPerUnitNorm = availableHeightNorm / totalVerticalUnits;
    fieldHeightNorm = heightPerUnitNorm;
    gapHeightNorm = heightPerUnitNorm * gapHeightFactor;

    % Check against minimum pixel height
    actualFieldHeightPixels = panelHeightPixels * fieldHeightNorm;
    if actualFieldHeightPixels < minFieldHeightPixels && panelHeightPixels > 0
        fieldHeightNorm = minFieldHeightPixels / panelHeightPixels; % Use minimum
        availableForGaps = availableHeightNorm - (numFields * fieldHeightNorm);
        if numGaps > 0 && availableForGaps > 0, gapHeightNorm = availableForGaps / numGaps; else gapHeightNorm = 0.005; end
        if (numFields * fieldHeightNorm + numGaps * gapHeightNorm) > availableHeightNorm, scaleFactor = availableHeightNorm / (numFields * fieldHeightNorm + numGaps * gapHeightNorm); fieldHeightNorm = fieldHeightNorm * scaleFactor; gapHeightNorm = gapHeightNorm * scaleFactor; end
    end

    % Position elements from top
    currentYNorm = 1.0 - innerTopMargin;
    for i = 1:length(fieldTags)
        tag = fieldTags{i};
        if strcmpi(tag, 'GAP'), currentYNorm = currentYNorm - gapHeightNorm; % Move down for gap
        else
            currentFieldBottomYNorm = currentYNorm - fieldHeightNorm; % Bottom Y of current field
            labelTag = ['label', tag]; valueTag = ['text', tag]; editTag = ['edit', tag];
            hLabel = findobj(hPanel, 'Tag', labelTag); hValue = findobj(hPanel, 'Tag', valueTag); hEdit = findobj(hPanel, 'Tag', editTag);
            if ~isempty(hLabel) && ishandle(hLabel), set(hLabel, 'Position', [innerMarginX, currentFieldBottomYNorm, labelW, fieldHeightNorm]); end
            hFieldToPos = []; if ~isempty(hValue) && ishandle(hValue), hFieldToPos = hValue; elseif ~isempty(hEdit) && ishandle(hEdit), hFieldToPos = hEdit; end
            if ~isempty(hFieldToPos), set(hFieldToPos, 'Position', [innerMarginX + labelW + gapX, currentFieldBottomYNorm, valueW, fieldHeightNorm]); end
            currentYNorm = currentFieldBottomYNorm; % Top Y for the next element
        end
    end
end

%% --- CORE CALLBACK FUNCTIONS ---

% =========================================================
% loadComponentsCallback: Carica dati componenti (CON DEBUG)
% =========================================================
function loadComponentsCallback(hObject, eventdata)
    hFigure = ancestor(hObject, 'figure'); handles = guidata(hFigure);
    localLogStatus(handles, 'Load Components button pressed.'); % Use local log
    set(hObject, 'Enable', 'off', 'String', 'Loading...'); drawnow;

    try
        % Get selections
        cartridgePopup = handles.popupmenuCartridge;
        bulletPopup = handles.popupmenuBullet;
        powderPopup = handles.popupmenuPowder;

        cartridgeList = get(cartridgePopup, 'String'); cartridgeIndex = get(cartridgePopup, 'Value');
        bulletList = get(bulletPopup, 'String'); bulletIndex = get(bulletPopup, 'Value');
        powderList = get(powderPopup, 'String'); powderIndex = get(powderPopup, 'Value');

        % --- DEBUGGING OUTPUT ---
        % (Keep debug prints if needed)
        % ...

        % Validate selections
        if ~iscell(cartridgeList) || ~iscell(bulletList) || ~iscell(powderList) || ...
           isempty(cartridgeList) || cartridgeIndex <= 1 || ...
           isempty(bulletList) || bulletIndex <= 1 || ...
           isempty(powderList) || powderIndex <= 1
             error('Select valid Cartridge, Bullet, and Powder from the dropdown menus.');
        end

        cartridgeName = cartridgeList{cartridgeIndex};
        bulletName = bulletList{bulletIndex};
        powderName = powderList{powderIndex};
        if contains(cartridgeName, {'(','No ','Err'}) || contains(bulletName, {'(','No ','Err'}) || contains(powderName, {'(','No ','Err'})
             error('Selected item name appears invalid.');
        end

        % --- Rest of the loading logic ---
        if exist('loadComponentData', 'file') ~= 2, error('loadComponentData function not found.'); end
        compData = loadComponentData(cartridgeName, bulletName, powderName);
        handles.componentData = compData;
        handles.isComponentsLoaded = true;
        handles.parametri = []; % Reset calculated parameters
        handles.isSimRun = false; % Reset simulation run flag
        handles.stressResults = []; % Reset stress results
        guidata(hFigure, handles);

        populateEditableFieldsFromData(hFigure); % Populate fields

        % Reset simulation parameters to defaults
        set(handles.editPropellantCharge_gr, 'String', '40.0');
        set(handles.editBarrelLength_in, 'String', '24.0');
        set(handles.editTwistRate_inPerTurn, 'String', '10.0');
        set(handles.editInitialTemperature_K, 'String', '293.15');
        set(handles.editAmbientPressure_Pa, 'String', '101325');
        set(handles.editShotStartPressure_MPa, 'String', '25.0');

        if exist('clearPlots','file')==2, clearPlots(hFigure);
        else, warning('Function clearPlots not found.'); end

        localLogStatus(handles, sprintf('Components loaded: %s/%s/%s. Select Barrel, verify/edit params & Apply.', cartridgeName, bulletName, powderName));

    catch ME
        handles.componentData = []; handles.isComponentsLoaded = false; handles.parametri = []; handles.isSimRun = false; handles.stressResults = [];
        guidata(hFigure, handles);
        localLogStatus(handles, ['ERROR loading components: ', ME.message]);
        uiwait(msgbox(sprintf('Error loading components:\n%s', ME.message), 'Load Error', 'error'));
        populateEditableFieldsFromData(hFigure); % Clear/reset fields
    end

    set(hObject, 'Enable', 'on', 'String', 'Load Selected Components'); % Re-enable button
    updateUIState(hFigure); % Update overall UI state
end


    % =========================================================
% runSimulationCallback: Executes simulation and updates barrel window
% =========================================================
function runSimulationCallback(hObject, ~)
% Executes the simulation and post-processing analysis.
% MODIFIED: Calls update function for barrel window if open.
% CORRECTED: Added missing call to displayStress2DWindow.

    hFigure = ancestor(hObject, 'figure');
    handles = guidata(hFigure);
    localLogStatus(handles, '"Run Simulation" pressed.');

    % --- Prerequisite Checks ---
    if ~isfield(handles,'isComponentsLoaded') || ~handles.isComponentsLoaded || ...
       ~isfield(handles,'barrelData') || isempty(handles.barrelData) || ...
       ~isfield(handles,'parametri') || isempty(handles.parametri)
        localLogStatus(handles, 'Error: Load components, select barrel, and apply parameters first.');
        uiwait(msgbox('Load components, select barrel, and press "Apply & Calculate Volumes" first.', 'Prerequisites Missing', 'warn'));
        return;
    end

    % --- Disable UI ---
    set(handles.pushbuttonLoadComponents, 'Enable', 'off');
    set(handles.pushbuttonSelectBarrel, 'Enable', 'off');
    set(handles.pushbuttonApplyAndCalc, 'Enable', 'off');
    set(hObject, 'Enable', 'off', 'String', 'Running...');
    if isfield(handles,'pushbuttonPowderParams') && ishandle(handles.pushbuttonPowderParams), set(handles.pushbuttonPowderParams,'Enable','off'); end
    drawnow;

    % --- Clear Plots and Reset State ---
    if exist('clearPlots','file')==2, clearPlots(hFigure); else warning('Function clearPlots not found.'); end
    handles.isSimRun = false; handles.risultati_sim = []; handles.stressResults = []; % Also clear stress results
    guidata(hFigure, handles);

    % --- Initialize variables for stress update ---
    stressResults = []; % Ensure it's initialized before the try block uses it
    barrelData = handles.barrelData; % Use barrel data already loaded in main GUI

    try
        localLogStatus(handles, 'Starting ODE solver...');
        if exist('runSimulation', 'file') ~= 2, error('runSimulation function not found.'); end
        if exist('simulationOdes', 'file') ~= 2, error('simulationOdes function not found.'); end

        % --- Run Simulation ---
        simResults = runSimulation(handles.parametri);

        % --- Process Results ---
        handles = guidata(hFigure); % Re-read handles in case they changed during sim pause
        handles.risultati_sim = simResults;

        if ~isempty(simResults) && isstruct(simResults) && isfield(simResults, 'timeS') && ~isempty(simResults.timeS) && length(simResults.timeS) > 1
            localLogStatus(handles, 'Simulation successful.');
            handles.isSimRun = true;
            guidata(hFigure, handles); % Save simulation results

            % --- Calculate Stresses ---
            if ~isempty(barrelData) && exist('calculateBarrelStresses', 'file') == 2
                try
                    localLogStatus(handles, 'Calculating barrel stresses...');
                    % Recalculate stressResults based on the latest simResults
                    stressResults = calculateBarrelStresses(handles.risultati_sim, barrelData); 
                    handles.stressResults = stressResults; % Store stress results in main handles
                    guidata(hFigure, handles); % Save stress results
                    localLogStatus(handles, sprintf('Stress calculated. Min SF: %.2f', stressResults.min_safety_factor));
                catch ME_stress
                    localLogStatus(handles, ['ERROR calculating stresses: ', ME_stress.message]);
                    handles.stressResults = []; % Clear on error
                    guidata(hFigure, handles);
                end
            elseif isempty(barrelData)
                 localLogStatus(handles, 'WARNING: Barrel data not loaded, skipping stress calculation.');
                 handles.stressResults = []; % Ensure it's clear if skipped
                 guidata(hFigure,handles);
            else
                 localLogStatus(handles, 'WARNING: calculateBarrelStresses function not found.');
                 handles.stressResults = []; % Ensure it's clear if skipped
                 guidata(hFigure,handles);
            end
            
            % --- NEED TO GET FRESH HANDLES AGAIN AFTER STORING STRESSRESULTS ---
            handles = guidata(hFigure); 
            % --- END FRESH HANDLES ---

            % --- Plotting and Basic Summary ---
            if exist('plotResults', 'file') == 2, plotResults(hFigure); % plotResults should now handle stressResults if present in handles
            else, localLogStatus(handles, 'WARN: plotResults function not found.'); end

            if exist('displaySummaryResults', 'file') == 2, displaySummaryResults(hFigure);
            else, localLogStatus(handles, 'WARN: displaySummaryResults function not found.'); end

            % --- Energy/BC Window ---
            try 
                % Check prerequisites for energy window again with fresh handles
                results_ok = isstruct(handles.risultati_sim) && all(isfield(handles.risultati_sim, {'timeS', 'projectileVelocityMps', 'angularVelocityRadps', 'frictionWorkJ', 'heatLossJ', 'gasTemperatureK', 'remainingPropellantMassKg', 'gasPressurePa'})) && ~isempty(handles.risultati_sim.timeS);
                % Note: BC/Form Factor removed from params_ok check in previous step
                params_ok = isstruct(handles.parametri) && all(isfield(handles.parametri, {'projMass_m', 'projMomentOfInertia_Ip', 'initialPropellantMass_m', 'specificGasConstant_R', 'specificHeatRatio_gamma', 'impetus_F', 'projWeight_gr', 'projDiameter_in'})); 
                barrel_ok = isstruct(handles.barrelData) && all(isfield(handles.barrelData, {'massKg', 'specificHeatJkgK'}));

                if results_ok && params_ok && barrel_ok
                    delta_T = NaN; balance_info = [];
                    if exist('calculateEnergyBalance', 'file') == 2
                        % Note: formFactor argument removed in previous step
                        [balance_info, delta_T] = calculateEnergyBalance(handles.risultati_sim, handles.parametri, handles.barrelData); 
                         % Note: BC calculation removed from calculateEnergyBalance
                    end
                     % Note: ballistic_coefficient_calculated argument removed in previous step
                    if exist('displayEnergyWindow', 'file') == 2, displayEnergyWindow(handles.risultati_sim, handles.parametri, delta_T); end 
                else 
                    localLogStatus(handles, 'Skipping energy window: Missing required data.'); 
                end
            catch ME_energyWin
                localLogStatus(handles, ['ERROR during energy balance/display: ', ME_energyWin.message]); 
                fprintf(2, 'Stack trace for energy/BC error:\n'); disp(ME_energyWin.getReport); 
            end

            % --- >>> BLOCCO AGGIUNTO/CORRETTO <<< ---
            % Visualizza la finestra separata degli stress 2D
            % Usa handles.stressResults aggiornato sopra
            if isfield(handles, 'stressResults') && ~isempty(handles.stressResults) && exist('displayStress2DWindow', 'file') == 2
                try
                    % Assicurati che anche barrelData sia disponibile (già verificato sopra)
                    if isfield(handles, 'barrelData') && ~isempty(handles.barrelData)
                        localLogStatus(handles, 'Opening 2D stress distribution window...'); % Messaggio log corretto
                        displayStress2DWindow(handles.stressResults, handles.barrelData); % Chiama la funzione corretta
                    else
                        % Questo caso non dovrebbe accadere se il calcolo stress è avvenuto
                        localLogStatus(handles, 'Skipping 2D stress window: Barrel data missing unexpectedly.');
                    end
                catch ME_dispStress
                    localLogStatus(handles, ['Error displaying 2D stress window: ', ME_dispStress.message]); % Messaggio log corretto
                end
            elseif ~isfield(handles, 'stressResults') || isempty(handles.stressResults)
                 % Questo log viene raggiunto se il calcolo stress è fallito o saltato
                 localLogStatus(handles, 'Skipping 2D stress window (no stress results calculated).'); 
            else % Function displayStress2DWindow non trovata
                 localLogStatus(handles, 'Skipping 2D stress window (displayStress2DWindow.m not found).');
            end
            % --- >>> FINE BLOCCO AGGIUNTO/CORRETTO <<< ---

            % --- >>> NEW: Update Barrel Window Display <<< ---
            if isfield(handles, 'barrelSelectionFigure') && ~isempty(handles.barrelSelectionFigure) && ishandle(handles.barrelSelectionFigure)
                localLogStatus(handles, 'Attempting to update barrel window with stress results...');
                hBarrelFig = handles.barrelSelectionFigure;
                barrelHandles = guidata(hBarrelFig);
                if ~isempty(barrelHandles) && isfield(barrelHandles, 'updateBarrelStressDisplay') && isa(barrelHandles.updateBarrelStressDisplay, 'function_handle')
                    if ~isempty(handles.stressResults) && ~isempty(handles.barrelData)
                        try
                            barrelHandles.updateBarrelStressDisplay(hBarrelFig, handles.stressResults, handles.barrelData);
                            localLogStatus(handles, 'Barrel window display updated.');
                        catch ME_update_barrel
                             localLogStatus(handles, ['ERROR calling updateBarrelStressDisplay: ', ME_update_barrel.message]);
                             fprintf(2,'Barrel Window Update Error Stack:\n'); disp(ME_update_barrel.getReport);
                        end
                    else
                         localLogStatus(handles, 'Skipping barrel window update (stress or barrel data missing).');
                    end
                else
                     localLogStatus(handles, 'WARNING: Could not find update function handle in barrel window guidata.');
                end
            else
                 localLogStatus(handles, 'Barrel window not open or handle invalid, skipping update.');
            end
            % --- >>> END NEW SECTION <<< ---

            % --- Save Results ---
             if exist('analyzeAndSaveResults', 'file') == 2
                 try 
                     analyzeAndSaveResults(handles.risultati_sim, handles.parametri); % Pass main handles structs
                 catch ME_save
                     localLogStatus(handles, ['ERROR saving results: ', ME_save.message]); 
                 end
             else 
                 localLogStatus(handles, 'WARN: analyzeAndSaveResults function not found.'); 
             end

        else % Simulation ran but produced no valid results
            localLogStatus(handles, 'Simulation ran but produced no valid results.');
            uiwait(msgbox('Simulation did not produce valid results.', 'Sim Error', 'warn'));
            handles.risultati_sim = []; handles.isSimRun = false; handles.stressResults = [];
            guidata(hFigure, handles);
        end

    catch ME % Error during runSimulation or analysis
        handles = guidata(hFigure); % Get fresh handles
        handles.risultati_sim = []; handles.isSimRun = false; handles.stressResults = [];
        guidata(hFigure, handles); % Save cleared state
        localLogStatus(handles, ['ERROR during simulation run: ', ME.message]);
        fprintf(2, 'Stack trace for simulation error:\n'); disp(ME.getReport);
        uiwait(msgbox(sprintf('Simulation Error:\n%s\nCheck log/console.', ME.message), 'Execution Error', 'error'));
    end

    % --- Re-enable UI ---
    try
        set(handles.pushbuttonLoadComponents, 'Enable', 'on');
        set(handles.pushbuttonSelectBarrel, 'Enable', 'on');
        componentsLoaded = isfield(handles,'isComponentsLoaded') && handles.isComponentsLoaded;
        barrelLoaded = isfield(handles,'barrelData') && ~isempty(handles.barrelData);
        set(handles.pushbuttonApplyAndCalc, 'Enable', iff(componentsLoaded && barrelLoaded,'on','off'));
        set(hObject, 'Enable', 'on', 'String', 'Run Simulation');
        if isfield(handles,'pushbuttonPowderParams') && ishandle(handles.pushbuttonPowderParams), set(handles.pushbuttonPowderParams,'Enable', iff(componentsLoaded,'on','off')); end
        updateUIState(hFigure);
    catch ME_ui_enable, localLogStatus(handles, ['Error re-enabling UI: ' ME_ui_enable.message]); end

end % End runSimulationCallback



function selectBarrelCallback(hObject, ~)
% Opens the separate barrel selection GUI window.
    hFigure = ancestor(hObject, 'figure'); handles = guidata(hFigure);
    localLogStatus(handles, 'Select Barrel button pressed.'); barrelFigField = 'barrelSelectionFigure';
    if isfield(handles, barrelFigField) && ~isempty(handles.(barrelFigField)) && ishandle(handles.(barrelFigField))
        localLogStatus(handles, 'Barrel selection window already exists. Bringing to front.'); figure(handles.(barrelFigField)); set(handles.(barrelFigField), 'Visible', 'on');
    else
            % Window doesn't exist or is invalid: Create it
            localLogStatus(handles, 'Creating new barrel selection window...'); 
            if exist('barrelSelectionGUI', 'file') == 2
                try
                    hBarrelSelFig = barrelSelectionGUI(hFigure); % Create the GUI [cite: 530]
                    
                    % --- AGGIUNTA: Controlla se la creazione ha avuto successo ---
                    if ishandle(hBarrelSelFig) 
                        handles.(barrelFigField) = hBarrelSelFig; % Store its handle [cite: 531]
                        guidata(hFigure, handles); % Save updated handles [cite: 532]
                        
                        % --- AGGIUNTA: Rendi la finestra visibile SUBITO ---
                        set(hBarrelSelFig, 'Visible', 'on'); 
                        figure(hBarrelSelFig); % Opzionale: porta anche in primo piano
                        % --- FINE AGGIUNTE ---
                    else
                         % Log o mostra errore se la creazione fallisce
                         localLogStatus(handles, 'ERROR: barrelSelectionGUI returned an invalid handle.');
                         uiwait(msgbox('Failed to create barrel selection window.', 'Error', 'error'));
                    end
                    % --- FINE AGGIUNTA ---

                catch ME_createBarrelGUI
                    localLogStatus(handles, ['ERROR creating barrel selection GUI: ', ME_createBarrelGUI.message]); [cite: 533]
                    uiwait(msgbox('Failed to create barrel selection window.', 'Error', 'error')); [cite: 533]
                end
            else
                 localLogStatus(handles, 'ERROR: barrelSelectionGUI.m function not found!'); [cite: 534]
                 uiwait(msgbox('Barrel selection GUI file is missing.', 'Error', 'error')); [cite: 534]
            end
    end
end

function openPowderParamsCallback_Wrapper(hObject, ~)
% Wrapper to call the external openPowderParamsCallback function.
    hFigure = ancestor(hObject, 'figure');
    if exist('openPowderParamsCallback', 'file') == 2
        try openPowderParamsCallback(hObject, []); catch ME_openPowder, localLogStatus(guidata(hFigure), ['ERROR calling openPowderParamsCallback: ', ME_openPowder.message]); fprintf(2, 'Stack trace:\n'); disp(ME_openPowder.getReport); uiwait(msgbox(sprintf('Error opening powder window:\n%s', ME_openPowder.message), 'Window Error', 'error')); end
    else localLogStatus(guidata(hFigure), 'ERROR: openPowderParamsCallback.m function not found!'); uiwait(msgbox('Powder parameters window function is missing.', 'Error', 'error')); end
end

%% --- LOCAL HELPER FUNCTIONS ---

function val = getfield_safe(s, field, default)
% Returns field value or default if field doesn't exist or is empty.
    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end

function out = iff(condition, trueVal, falseVal)
% Simple inline-if function substitute.
    if condition, out = trueVal; else, out = falseVal; end
end

function localLogStatus(handles, message)
% Logs message to the log window if available, otherwise to console.
% Uses handles struct passed from the calling function.
    if isfield(handles, 'logEditBox') && ~isempty(handles.logEditBox) && ishandle(handles.logEditBox)
        try
            hLog = handles.logEditBox; hLogFig = ancestor(hLog, 'figure'); timestamp = datestr(now, 'HH:MM:SS'); currentLog = get(hLog, 'String'); if isempty(currentLog), currentLog = {}; elseif ~iscell(currentLog), currentLog = {currentLog}; end
            if ischar(message), messageLines = strsplit(message, '\n'); elseif iscell(message), messageLines = message; else messageLines = {'Log Err: Invalid msg type'}; end
            newLogEntries = cell(sum(~cellfun('isempty', strtrim(messageLines))), 1); count = 0; for i = 1:length(messageLines), lineTxt = strtrim(messageLines{i}); if ~isempty(lineTxt), count=count+1; newLogEntries{count} = sprintf('[%s] %s', timestamp, lineTxt); end; end
            updatedLog = [currentLog; newLogEntries(1:count)]; maxLines = 200; if length(updatedLog) > maxLines, updatedLog = updatedLog(end-maxLines+1:end); end
            set(hLog, 'String', updatedLog); set(hLog, 'Value', length(updatedLog)); if ishandle(hLogFig) && strcmp(get(hLogFig,'Visible'), 'off'), set(hLogFig,'Visible','on'); end; drawnow limitrate;
        catch ME_log, fprintf(2, 'Internal localLogStatus error: %s\n', ME_log.message); disp(['LOG (Error): ', message]); end
    else disp(['LOG (No Win): ', message]); end
end

% =========================================================
% Wrapper per chiamare applyAndCalcCallback esterno
% =========================================================
function applyAndCalcCallback_Wrapper(hObject, eventdata)
    hFigure = ancestor(hObject, 'figure');
    handles = guidata(hFigure); % Get handles just in case
    localLogStatus(handles, 'Apply & Calc Wrapper called...'); % Use local log
    set(hObject, 'Enable', 'off', 'String', 'Calculating...'); drawnow; % Disable button

    calculation_successful = false; % Flag
    try
        if exist('applyAndCalcCallback', 'file') == 2
             calculation_successful = applyAndCalcCallback(hFigure); % Chiama funzione esterna
        else error('Function applyAndCalcCallback.m not found in path.'); end

        if calculation_successful
             localLogStatus(handles, 'Apply & Calc successful.');
             % ADDED: Store selectedBarrelName in parameters struct
             handles = guidata(hFigure); % Get updated handles (containing parameters)
             if isfield(handles, 'selectedBarrelName') && ~isempty(handles.selectedBarrelName)
                 if isfield(handles, 'parametri') && ~isempty(handles.parametri)
                     handles.parametri.selectedBarrelName = handles.selectedBarrelName;
                     guidata(hFigure, handles); % Save updated parameters
                     localLogStatus(handles, sprintf('Stored selected barrel name "%s" in parameters.', handles.selectedBarrelName));
                 else
                     localLogStatus(handles, 'WARNING: Parameters struct not found after successful Apply&Calc. Cannot store barrel name.');
                 end
             else
                  localLogStatus(handles, 'WARNING: No selected barrel name found in handles after Apply&Calc.');
             end
        else
            localLogStatus(handles, 'Apply & Calc failed (see previous errors).');
            handles = guidata(hFigure); handles.parametri = []; handles.stressResults = []; guidata(hFigure, handles); % Clear params and stress results
        end
    catch ME
        localLogStatus(handles, ['CRITICAL ERROR during external Apply & Calc call: ', ME.message]);
        fprintf(2, 'Stack trace:\n'); disp(ME.getReport);
        uiwait(msgbox(sprintf('Critical Error:\n%s', ME.message), 'Calculation Error', 'error'));
        try handles_err = guidata(hFigure); handles_err.parametri = []; handles_err.stressResults = []; guidata(hFigure, handles_err); catch, end; calculation_successful = false;
    end

    % Riabilita il bottone e aggiorna UI state
    try
        set(hObject, 'Enable', 'on', 'String', 'Apply & Calculate Volumes');
        if exist('updateUIState','file') == 2, updateUIState(hFigure); else warning('updateUIState function not found.'); end
    catch ME_reenable, localLogStatus(handles, ['Error re-enabling Apply button: ' ME_reenable.message]); end
end % Fine applyAndCalcCallback_Wrapper


% --- populateEditableFieldsFromData ---
% (Keep implementation from previous response)
function populateEditableFieldsFromData(hFigure)
    handles = guidata(hFigure); if ~isfield(handles, 'componentData') || isempty(handles.componentData), return; end; cd = handles.componentData; C = handles.CONV; localLogStatus(handles, 'Populating editable fields from loaded data...');
    try
        set(handles.editBulletName, 'String', getfield_safe(cd.bullet, 'bulletName', 'N/A')); mass_gr = getfield_safe(cd.bullet, 'mass_kg', NaN) * C.GRAINS_PER_KG; set(handles.editBulletMass_gr, 'String', iff(isnan(mass_gr), '', sprintf('%.1f', mass_gr))); length_in = getfield_safe(cd.bullet, 'length_m', NaN) * C.INCHES_PER_METER; set(handles.editBulletLength_in, 'String', iff(isnan(length_in), '', sprintf('%.3f', length_in))); diam_in = getfield_safe(cd.bullet, 'diameter_m', NaN) * C.INCHES_PER_METER; set(handles.editBulletDiameter_in, 'String', iff(isnan(diam_in), '', sprintf('%.3f', diam_in)));
        set(handles.editCartridgeName, 'String', getfield_safe(cd.cartridge, 'cartridgeName', 'N/A')); caseLen_in = getfield_safe(cd.cartridge, 'caseLength_m', NaN) * C.INCHES_PER_METER; set(handles.editCaseLength_in, 'String', iff(isnan(caseLen_in), '', sprintf('%.3f', caseLen_in))); maxCap_grH2O = getfield_safe(cd.cartridge, 'maxCaseCapacity_m3', NaN) * C.GRH2O_PER_M3; set(handles.editMaxCaseCapacity_grH2O, 'String', iff(isnan(maxCap_grH2O), '', sprintf('%.2f', maxCap_grH2O))); boreDiam_in = getfield_safe(cd.cartridge, 'boreDiameter_m', NaN) * C.INCHES_PER_METER; set(handles.editBoreDiameter_in, 'String', iff(isnan(boreDiam_in), '', sprintf('%.3f', boreDiam_in))); seatingDepth_m = getfield_safe(cd.bullet, 'seatingDepth_m', NaN); seatingDepth_in = seatingDepth_m * C.INCHES_PER_METER; set(handles.editSeatingDepth_in, 'String', iff(isnan(seatingDepth_in), '', sprintf('%.3f', seatingDepth_in)));
        set(handles.textCartridgeOAL_in, 'String', '---'); set(handles.textVolumeOccupiedBullet_grH2O, 'String', '---'); set(handles.textQLUseableCaseCapacity_grH2O, 'String', '---'); set(handles.textInitialFreeVolumeV0_grH2O, 'String', '---'); set(handles.textBulletTravel_in, 'String', '---');
    catch ME_populateEdit, localLogStatus(handles, ['ERROR populating editable fields: ', ME_populateEdit.message]); end
end

% --- updateUIState (ensure this function exists and is accessible) ---
% Placeholder if external file is missing
function updateUIState(hFigure)
    disp('updateUIState called (placeholder - ensure external file exists)');
    % Add basic logic here if needed as fallback
    handles = guidata(hFigure);
    componentsLoaded = isfield(handles, 'isComponentsLoaded') && handles.isComponentsLoaded;
    barrelLoaded = isfield(handles, 'barrelData') && ~isempty(handles.barrelData);
    paramsCalculated = isfield(handles, 'parametri') && ~isempty(handles.parametri);

    set(handles.pushbuttonApplyAndCalc, 'Enable', iff(componentsLoaded && barrelLoaded, 'on', 'off'));
    set(handles.pushbuttonRunSimulation, 'Enable', iff(paramsCalculated, 'on', 'off'));
    % Add more enable/disable logic as needed
end

