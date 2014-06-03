function [value,isterminal,direction] = bgcEventThrustDown(t,y,zdrop,varargin)

% Decompose state vector
zt = y(1);
z = y(2);
 
value = z-zdrop; % increases during descent.
isterminal = true;
direction = +1;