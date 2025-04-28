% gui/createLogWindow.m
function hLogEditBox = createLogWindow(mainFigPos)
% Crea una finestra separata per visualizzare i messaggi di log.
% INPUT:
%   mainFigPos: Posizione normalizzata della finestra principale [left bottom width height]
%               usata per posizionare la finestra di log accanto.
% OUTPUT:
%   hLogEditBox: Handle alla casella di testo editabile nella finestra di log.

logFigWidth = 0.3; % Larghezza normalizzata rispetto allo schermo
logFigHeight = 0.4; % Altezza normalizzata rispetto allo schermo

% Calcola posizione per metterla a destra della finestra principale
logFigLeft = mainFigPos(1) + mainFigPos(3) + 0.01;
logFigBottom = mainFigPos(2) + mainFigPos(4) - logFigHeight;
% Assicura che rimanga nello schermo
screenSize = get(0, 'ScreenSize'); % [1 1 width height] in pixels
screenWidthNorm = screenSize(3); screenHeightNorm = screenSize(4); % These are not normalized, fix needed if using normalized screen units
% Using normalized figure units requires careful handling relative to screen size
% For simplicity, let's position relative to main window in normalized units
if logFigLeft + logFigWidth > 0.98 % Evita di uscire a destra
    logFigLeft = mainFigPos(1) - logFigWidth - 0.01;
end
if logFigLeft < 0.01, logFigLeft = 0.01; end % Evita di uscire a sinistra
if logFigBottom < 0.05, logFigBottom = 0.05; end % Evita di uscire in basso

logPosition = [logFigLeft, logFigBottom, logFigWidth, logFigHeight];

hFig = figure('Name', 'Simulator Log', ...
              'Units', 'normalized', ...
              'Position', logPosition, ...
              'NumberTitle', 'off', ...
              'MenuBar', 'none', ...
              'ToolBar', 'none', ...
              'Resize', 'on', ...
              'CloseRequestFcn', @logWindowCloseReq, ... % Impedisce chiusura diretta
              'Tag', 'logWindowFigure', ... % Tag per trovarla
              'Visible', 'on');

hLogEditBox = uicontrol('Parent', hFig, ...
                        'Style', 'edit', ...
                        'Units', 'normalized', ...
                        'Position', [0.02 0.02 0.96 0.96], ... % Quasi tutta la finestra
                        'Min', 0, 'Max', 2, ... % Multi-line
                        'Enable', 'inactive', ... % Read-only look
                        'HorizontalAlignment', 'left', ...
                        'BackgroundColor', [1 1 1], ... % Sfondo bianco
                        'FontName', 'Monospaced', ... % Font leggibile per log
                        'String', {'--- Log Start ---'}, ...
                        'Tag', 'logEditBox');

% Nascondi la finestra invece di chiuderla
    function logWindowCloseReq(src, ~)
        set(src, 'Visible', 'off');
        disp('Log window hidden.');
    end

end