function [value,isterminal,direction] = bgcEventControl(t,y,zFilter,zTol,tControl,tLockout,varargin)
% For use with ODE45.
%
% Revision History
% 2012-12-30    mvj    Created.
% 2013-01-02    mvj    Added time-based disengagement behavior.  Kludgey.


% Use timeout to disengage controller.
persistent tEngage;

% Decompose state vector
zt = y(1);
z = y(2);
% Engage controller on upcast when within depth bands of interest.
if zt < 0
  zValue = min(abs(z-zFilter)-zTol);  % < 0 when within a control band.
else
  zValue = 1;
end

[~,ii] = min(abs(z-zFilter));
[zt z zFilter(ii) zValue]

% Mark engagement time.
if zValue < 0 && isempty(tEngage)
  tEngage = t;
elseif (t-tEngage) - tControl - tLockout > 0
  tEngage = [];
end

% Engagement of controller is based on z criteria.  
% Disengage after a period of activation.
if isempty(tEngage)
  value = zValue;
else
  value = (t-tEngage) - tControl;
end
  
isterminal = true;
direction = 0;
