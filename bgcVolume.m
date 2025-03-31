function V = bgcVolume(Vo,alpha,chi,dtheta,dp,varargin)
%
% Linearized equation of state for a solid or liquid.
%
% Vo: initial volume.
% alpha: linear coefficient of thermal expansion
% chi: volumetric bulk compressibility
% dtheta: change from initial temperature
% dp: change from initial pressure.
%
% Revision History
% 2013-01-03    mvj    Created.
% 2025-03-31    mvj    make compatible with other equations of state requiring more parameters.

% Compute volume at this temperature and pressure
V = Vo*(1 + 3*alpha*dtheta - chi*dp);
