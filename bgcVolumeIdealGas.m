function V = bgcVolumeIdealGas(Vo,alpha,chi,dtheta,dP,thetao,Po,varargin)
%
% Vo, To, Po: initial conditions.
% T, P: in situ temperature (K) and pressure (Pa)
%
% Applies the ideal gas law: PV = nRT
%
% alpha and chi are ignored. (Retained for backwards compatibility.)
%
% Revision History
% 2025-03-31    mvj    Created.

% Compute volume at this temperature and pressure
theta = dtheta + thetao;
P = dP + Po;
V = theta./P .* Vo.*Po./thetao;

% Optionally handle a lockout volume.
if nargin > 8
    izg = varargin{1};
    VLockouts = varargin{2}{2};
    izg = min(izg,length(VLockouts));  % hack
    V = max(VLockouts(izg),V);
elseif nargin > 7
    VLockout = varargin{1}{2};
    V = max(VLockout,V);
end

