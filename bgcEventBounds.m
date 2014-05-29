function [value,isterminal,direction] = bgcEventBounds(t,y,zSurface,zSeafloor)
% For use with ODE45.
%
% Revision History
% 2013-01-02    mvj    Created.

% Decompose state vector
zt = y(1);
z = y(2);

% Sea surface and seafloor
value = min(zSeafloor - z,z - zSurface);

isterminal = true;
direction = 0;
