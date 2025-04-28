% =========================================================
% Internal Ballistics Simulator - Refactored Version
% =========================================================
% core/barrelExitEvent.m
% Event function for the ODE solver (ode45, ode15s, etc.).
% Detects when the projectile reaches the end of the barrel.
% =========================================================

function [value, isterminal, direction] = barrelExitEvent(t, y, barrelLength)
% INPUTS:
%   t: Current time (scalar) - provided by ODE solver
%   y: Current state vector (column) - provided by ODE solver
%      y = [remainingPropellantMass; gasTemperature; projectilePosition;
%           projectileVelocity; angularVelocity; frictionWork; heatLoss]
%   barrelLength: Length of the barrel [m] - passed as an additional parameter from odeset
% OUTPUTS:
%   value: Value of the event function. Event occurs when value = 0.
%          Here: value = projectilePosition - barrelLength
%   isterminal: 1 if integration should terminate at the event, 0 otherwise.
%   direction: 0 detects zeros regardless of direction (default)
%              1 detects only increasing zeros (value crosses from - to +)
%             -1 detects only decreasing zeros (value crosses from + to -)

    % Index of the projectile position (projectilePosition) in the state vector y
    PROJECTILE_POSITION_INDEX = 3; % Ensure this index matches the state vector definition

    % Extract the current projectile position
    currentProjectilePosition = y(PROJECTILE_POSITION_INDEX);

    % Calculate the event function value
    % The event occurs when projectilePosition reaches barrelLength.
    value = currentProjectilePosition - barrelLength;

    % Terminate the integration when the event occurs
    isterminal = 1;

    % Detect the event only when the position is increasing (crossing barrelLength)
    direction = 1;

end % End of function barrelExitEvent