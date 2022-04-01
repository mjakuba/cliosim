function omega = bgcBruntVaisala(prm)
% Undamped frequency of oscillation for profiler at neutral depth.  Assumes a 
% stable profiler.  Output is in rad/s.
%
% Revision History
% 2012-12-28    mvj    Created.

% Find equilibrium depth.  First step finds interval within which equilibrium depth lies.
rhof = (prm.m/prm.V)./(1+prm.alpha*(prm.profile.theta-prm.theta)-prm.chi*(prm.profile.p-prm.const.atm));
izabove = find(rhof >= prm.profile.rho,1,'last');
rhof
izabove
zabove = prm.profile.z(izabove);
zbelow = prm.profile.z(izabove+1);
zabove
zbelow
[rhoabove,thetaabove,pabove,drhodz,dthetadz,dpdz] = bgcInSitu(zabove,prm.profile);

% Second step uses linear profile within interval to compute the equilibrium
% depth analytically.
% @@@ a is the problem.  this should still work for alpha, chi = 0
a = drhodz*dthetadz*prm.alpha - drhodz*dpdz*prm.chi;
b = drhodz;
c = rhoabove*(1+prm.alpha*(thetaabove-prm.theta) - prm.chi*(pabove-prm.const.atm)) - prm.m/prm.V;
z = zabove + (-b + [1 -1]*sqrt(b^2-4*a*c))/2/a;
a
b
c
z
z = z(z >= zabove & z < zbelow);

% Compute in situ properties.
[rho,theta,p,drhodz,dthetadz,dpdz] = bgcInSitu(z,prm.profile);

% Undamped freqency of oscillation (linearized).  See notes, 2012-Dec-28
omega = sqrt(prm.const.g/rho*(drhodz - prm.chi*rho^2*prm.const.g + prm.alpha*dthetadz*rho));

