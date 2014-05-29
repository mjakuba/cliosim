function V = bgcVolume(Vo,alpha,chi,dtheta,dp)
%
% Vo: initial volume.
% alpha: linear coefficient of thermal expansion
% chi: volumetric bulk compressibility
% dtheta: change from initial temperature
% dp: change from initial pressure.
%
% Revision History
% 2013-01-03    mvj    Created.

% Compute volume at this temperature and pressure
V = Vo*(1 + 3*alpha*dtheta - chi*dp);
