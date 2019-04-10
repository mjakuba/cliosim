function [value,isterminal,direction] = bgcEventFilter(t,y,zFilter,zTol,ztTol,tPump,tLockout)
% For use with ODE45.
%
% Revision History
% 2012-12-30    mvj    Created.
% 2013-01-02    mvj    Added timeout-based disengagement.
% 2019-04-10    mvj    This is inherently brittle.  It attempts to guess the mission plan excuted in bgcF based on state.


% Use timeout to disengage pumping.
persistent tEngage

% Decompose state vector
zt = y(1);
z = y(2);

% Filter in specified depth bands on upcast. 
% Supports vector zFilter for multiple filtering depths.
if abs(zt) < ztTol
  zValue = min(abs(z-zFilter)-zTol);
else
  zValue = 1;
end

% Mark engagement time.
if zValue < 0 && isempty(tEngage)
  tEngage = t;
elseif (t-tEngage) - tPump - tLockout > 0  % prevent immediate re-triggering with lockout.
  tEngage = [];
end

% Engagement of pump is based on z,zt criteria.  
% Disengage after a period of activation.
if isempty(tEngage)
  value = zValue;
else
  value = (t-tEngage) - tPump;
end

isterminal = true;
direction = 0;
