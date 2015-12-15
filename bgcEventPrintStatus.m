function [value,isterminal,direction] = bgcEventPrintStatus(t,y,DT)
% For use with ODE45.
%
% Revision History
% 2015-12-15    mvj    Created.

% Decompose state vector
zt = y(1);
z = y(2);

fprintf(1,'t = %0.1f s z = %0.1f zt = %0.3f\n',t,z,zt);

value = 1.0;
isterminal = false;
direction = 0;
