function [value,isterminal,direction] = bgcEventNone(t,y)
% For use with ODE45.  An event that will never be true.
%
% Revision History
% 2012-12-28    mvj    Created.

value = 1;
isterminal = 1;
direction = 0;

