function prm = bgcConst()
% Default constants.
%
%
% Revision History
% 2014-09-08    mvj    Moved from vehicle/mission definition files.
% 2022-03-12    mvj    Add default latitude for pressure/depth conversion.


prm.const.g = 9.81; % [m/s^2]
prm.const.nu = 1.83e-6; % [m^2/s, 0 C]
prm.const.mu = 1.88e-3; % [Pa s, 0 C]
prm.const.atm = 101325; % [Pa] atmospheric pressure
prm.const.lat = 0;
