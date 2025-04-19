function [value,isterminal,direction] = bgcEventBounds(t,y,zSurface,zSeafloor,varargin)
% For use with ODE45.
%
% Revision History
% 2013-01-02    mvj    Created.
% 2025-04-16    mvj    Add option to continue sim---logic is in bgc.m; added varargin to pass additional parameters.


% Decompose state vector
zt = y(1);
z = y(2);

% Sea surface and seafloor
value = min(zSeafloor - z,z - zSurface);

isterminal = true;
direction = 0;
