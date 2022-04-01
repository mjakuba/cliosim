function [rho,theta,p,drhodz,dthetadz,dpdz,S,dSdz] = bgcInSitu(z,profile)
% Compute basic in situ properties at depth z for a given profile.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2012-12-28    mvj    Added output parameters.
% 2012-12-28    mvj    defined slopes at a profile depth z as the slope immediately below z.
% 2022-03-04    mvj    Add in situ salinity as an output.


izabove = find(z >= profile.z,1,'last');
izbelow = find(z < profile.z(izabove:end),1)+izabove-1;
assert(~isempty(izbelow.*izabove),sprintf('Profiler is outside of profile data (z = %.1f)',z));
iz = [izabove izbelow];
dthetadz = diff(profile.theta(iz))/diff(profile.z(iz));
drhodz = diff(profile.rho(iz))/diff(profile.z(iz));
dpdz = diff(profile.p(iz))/diff(profile.z(iz));
rho = profile.rho(izabove) + drhodz*(z-profile.z(izabove));
theta = profile.theta(izabove) + dthetadz*(z-profile.z(izabove));
p = profile.p(izabove) + dpdz*(z-profile.z(izabove));

dSdz = diff(profile.S(iz))/diff(profile.z(iz));
S = profile.S(izabove) + dSdz*(z-profile.z(izabove));
