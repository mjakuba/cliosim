function Vfo = bgcNeutral(z,profile,prmc,prmf)
% Vfo = bgcNeutral(z,profile,prmc,prmf)  Compute the volume V of flotation with 
% parameters prmf to make the profiler neutral at depth z.  prmc  denotes bulk
% parameters for the remaining components of the profiler.  bgcParam should be
% used to create the prm* structures.
%
% Revision History
% 2012-12-27    mvj    Created.

% Compute in situ properties.
[rho,theta,p] = bgcInSitu(z,profile);

% Compute displacement of the components at in situ temperature and pressure
Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prmc.theta),(p-prmc.const.atm));

% Compute the required buoyancy.
Zc = prmc.const.g*prmc.m - Vc*rho*prmc.const.g;
assert(Zc > 0,'Float is buoyant without any flotation!');
  
% Compute the density of the floatation at in situ temperature and pressure
% Nominal temperature and pressure are assumed to be the same as for other components.
vf = bgcVolume(1,prmf.alpha,prmf.chi,(theta-prmc.theta),(p-prmc.const.atm));
rhof = prmf.rho/vf;

% Compute the volume of flotation required at pressure.
Vf = -Zc/prmc.const.g/(rhof - rho);

% Compute required volume at surface.
Vfo = Vf/vf;

% Cross-check.  Net buoyancy at depth should be numerically zero.
alpha = (prmc.alpha*prmc.V + Vfo*prmf.alpha)/(prmc.V+Vfo);
chi = (prmc.chi*prmc.V + Vfo*prmf.chi)/(prmc.V+Vfo);
Vo = prmc.V+Vfo;
V = bgcVolume(Vo,alpha,chi,(theta-prmc.theta),(p-prmc.const.atm));
m = prmc.m + Vfo*prmf.rho;
Z = prmc.const.g*m - V*rho*prmc.const.g;
assert(abs(Z)<1e-9,'Float is not neutral at desired depth!');

