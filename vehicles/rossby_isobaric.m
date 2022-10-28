function [prm,components] = bgcParam
% isobaric Rossby float, idealized
%
% Revision History
% 2022-03-02    mvj    Created.

% o use profile from Heather for off Madagascar, but first use an isothermal and isohaline profile
% o use some kind of controller to get 

% Mission parameters
BALLAST_DEPTH = 100;

% Load conversion constants.
bgcConversions;

% Still need this file for things that aren't in the model - fluids etc.
bgcMatl;

% Initialize parameter structure with constants, override as necessary.
prm = bgcConst();

% Initial conditions and other things related to the solver.
prm.solver.zo = BALLAST_DEPTH; % [m]
prm.solver.zto = -0.5; % [m/s]
prm.solver.thetao = NaN;  % [K] not implemented.
prm.solver.tend = 3600*10; % s

% The background water column profile.
%prm.profile = bgcProfile(prm.const,'isopycnal/isothermal');
%prm.profile =  bgcProfile(prm.const,'pycnoclinic/isohaline');
prm.profile =  bgcProfile(prm.const,'pycnoclinic/positive-upward haloclinic');
%prm.profile =  bgcProfile(prm.const,'pycnoclinic/positive-downward haloclinic');
%prm.profile = bgcProfile(prm.const,'levitus82World');
% ARGO profiles from Agulhas Current.  
%argo = load('~/Dropbox/jakuba/veh/rafos_compressee/profiles/ArgoProfileMatrices_forMJ.mat');
%argo.ii = 56; % Profiles 56 and 83 decently bracket the variability.
%prm.profile = bgcProfile(prm.const,argo.argoPS(:,argo.ii),argo.argoTE(:,argo.ii),argo.argoPR(:,argo.ii));

% Move the water column.  The simulation does not accurately simulated moving water, so the results
% really only apply to the steady state.
prm.profile.displacement.type = 'step'
prm.profile.displacement.step.tstep = 2*3600;
prm.profile.displacement.step.step = 50;

% The BGC profiler as a whole has:
% * a temperature on deck, theta [K]
% * added mass, h [kg]
% * stokes drag coefficient, CDs [-]
% * skin friction coefficient, CDf [-], referenced to surface area
% * downcast quadratic drag coefficient, CDd [-], referenced to frontal area
% * upcast quadratic drag coefficient, CDu [-], referenced to frontal area
% * surface area, As [m^2]
% * frontal area, Af [m^2]
% * characteristic diameter, D [m]
% @@@ will need some damping for this to work at all.
prm.theta = 300; % [K]
prm.h = 10; % [kg]
prm.CDs = 0; % [-]
prm.CDf = 0.001; % [-]
prm.CDd = 0.1; % [-]  
prm.CDu = 0.1; % [-]
prm.As = 1.0; % [m^2]
prm.Af = 0.83; % [m^2]
prm.D = 0.5; % [m]

% Each component has:
% * a name 
% * mass, m [kg]
% * displacement, V [m^3], at 1 atm (zero if internal)
% * coefficient of thermal expansion, alpha [1/K]
% * an adiabatic compressibility, chi [1/Pa]
% * specific heat at constant pressur, cp [J/kg/K] (not used presently)
%

c = bgcInitComponent('Lumped vehicle');
c.rho = 1030;
c.V = pi*0.1^2*1;
c.m = c.V*c.rho;
%c.chi = 0;  % isobaric (really this means iso-in-situ-density)
%c.alpha = 0; % isothermal (optimal)
c.alpha = 0.6e-4; %0.5e-4;  % this is the linear alpha!  this gets slightly decreased by alpha=0 ballast added below.
%c.chi = 1/seawater.bulkModulus; % isopycnal to decent approximation.
[~,T_K,p_Pa,~,~,~,S_PSU,~]=bgcInSitu(BALLAST_DEPTH,prm.profile);
c.chi = 1/sw_seck(S_PSU,T_K-273.15,p_Pa*1e-4)/1e5; % isopycnal at target depth
c.chi = 0.98*c.chi;  % margin
c.cp = NaN;
prm.components = bgcAddComponent(c);

% compressee
% c = bgcInitComponent('Silicone Oil');
% c.rho = siliconeoil.density;
% c.V = 200/1000;  % roughly 200 l gives a net bulk modulus of 2.2 GPa (water) for k = 3.6e9.
% c.m = c.V*c.rho;
% c.alpha = siliconeoil.coeffThermalExpansion; % no info
% c.chi = 1/siliconeoil.bulkModulus; % crude - this varies a lot.
% c.cp = siliconeoil.specificHeat;
% prm.components = bgcAddComponent(c,prm.components);


% ballasting.  This should be small relative to lumped vehicle above. 
% Add ballast/flotation as necessary to make vehicle neutral at 3000 m.
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
[rho,theta,p] = bgcInSitu(BALLAST_DEPTH,prm.profile);
Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); % Volume of the float at depth.
Zc = prm.const.g*prmc.m - Vc*rho*prm.const.g; % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
if Zc < 0
  fprintf(1,'Vehicle is positive.  Approx. %.1f kg margin ballast yields neutral at %.1f m\n',-Zc/prm.const.g,BALLAST_DEPTH);
else
  fprintf(1,'Vehicle is negative at %.1f m.  %.1f N additional floatation required.\n',BALLAST_DEPTH,Zc);
end

if Zc < 0
  f = bgcInitComponent('ballast (negatively buoyant)');
  f.rho = 10000;
  f.alpha = 0; 
  f.chi = 0; 
  margin = 1;
else
  f = bgcInitComponent('ballast (positively buoyant)');
  f.rho = 0;
  f.alpha = 0;
  f.chi = 0;
  margin = 1;  % this is the margin on the ballast, not on the whole vehicle.
end
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
f.V = bgcNeutral(BALLAST_DEPTH,prm.profile,prmc,f)*margin;
f.m = f.rho*f.V;
prm.components = bgcAddComponent(f,prm.components);

% Descent controller not used for floats but needs to appear in component list.
c = bgcInitComponent('descentController');
c.active = 1;  % this is actually not-active for this component only.  See bgcF.m
prm.components = bgcAddComponent(c,prm.components);

% Drop weight not used for floats but needs to appear in component list.
c = bgcInitComponent('drop weight');
c.active = 0; 
prm.components = bgcAddComponent(c,prm.components);

% Surface and seafloor.  For convenience these are massless components
c = bgcInitComponent('bounds');
c.eventf = @bgcEventBounds;
c.event_prm = {0,1000.0};
prm.components = bgcAddComponent(c,prm.components);

% Controller is always last component.  Not active for float sims but needs to appear in component list.
c = bgcInitComponent('controller');
c.active = 0;;
prm.components = bgcAddComponent(c,prm.components);


