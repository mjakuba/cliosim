function [value,isterminal,direction] = bgcEventDropWeight(t,y,zdrop)
% For use with ODE45.
%
% Revision History
% 2012-12-28    mvj    Created.

% Decompose state vector
zt = y(1);
z = y(2);

value = z-zdrop; % increases during descent.
isterminal = true;
direction = +1;
