function Vfo = bgcNeutral(z,profile,prmc,prmf)
% Vfo = bgcNeutral(z,profile,prmc,prmf)  Compute the volume V of flotation with 
% parameters prmf to make the profiler neutral at depth z.  prmc denotes effective
% lumped parameters for a portion of a vehicle, generated by passing all or part
% of a complete parameters structure (as generated by bgcParam) to bgcBulkParam:
%
% >> prmc = prm;
% >> [prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
% 
% prmf should be a single component substructure of a full simulation parameters
% definition, and generally should be the flotation, i.e. prm.components(n) 
% where prm.components(n).name = 'floatation'.
% 
%
% Revision History
% 2012-12-27    mvj    Created.
% 2014-10-08    mvj    Clarified usage.
% 2015-09-03    mvj    Altered to allow adding margin ballast if material is 
%                      negatively buoyant.


% Compute in situ properties.
[rho,theta,p] = bgcInSitu(z,profile);

% Compute displacement of the components at in situ temperature and pressure
Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prmc.theta),(p-prmc.const.atm));

% Compute the density of the floatation at in situ temperature and pressure
% Nominal temperature and pressure are assumed to be the same as for other components.
vf = bgcVolume(1,prmf.alpha,prmf.chi,(theta-prmc.theta),(p-prmc.const.atm));
rhof = prmf.rho/vf;

% Compute the required buoyancy.
Zc = prmc.const.g*prmc.m - Vc*rho*prmc.const.g;

if Zc < 0  % (N) buoyancy (<0 indicates float is positive, >0 float is negative).
  assert(rhof > rho,'Float is buoyant without any flotation!');  % negative ballast needed.
else
  assert(rhof < rho,'Float is negative without any ballast!');  % positive flotation needed.
end
  
% Compute the volume of flotation/ballast required at pressure.
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

