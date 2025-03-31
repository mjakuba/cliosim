function [prm] = bgcParam
% Armada SBIR phase 1 float concept.
%
% Focus is on mesoscale.
%
% Revision History
% 2025-03-31    mvj    Created.

% General approach:
% * coarsely model a TR-sized float
% * add a controller and passive buoyancy of appropriate size to effect near-surface buoyancy control
% * show that for a float less compressible than seawater it is possible to passively float at some neutral depth.
% * examine sensitivity to certain factors as follows:
% ** accuracy of pre-ballasting - how accurately can a particular depth be chosen?
% ** variation in lumped float parameters - use experimentally observed variation in RAFOS float parameters.
% ** variation in background profile - how isobaric will floats be?
% ** material choices

% TODO
% o get a representative ocean profile with mesoscale features.
% Shelf Research Fleet: https://scienceweb.whoi.edu/seasoar/cfrfwhoi/index_profiles.html

% Mission parameters.  Adapt Clio mission to float mission.  Relevant logic is in bgcF.m
% Ascent is triggered after all samples are complete (bgcF.m).
BALLAST_DEPTH = 100; 
CUTOFF_DEPTH = 5; % [m] initial descent to below unstable near-surface neutral depth
sampleDepths = BALLAST_DEPTH; % [m] Desired stable neutral depth.
sampleDepthTol = 5; % +/- [m] Sample timer starts once within this band.
sampleDepthRateTol = 0.1; % +/- [m/s] and depth rate below this figure.
sampleTime = 600; % [s] time to remain at sample depth.
sampleTimeLockout = 60; % [s] ?


% Load conversion constants.
bgcConversions;

% Still need this file for things that aren't in the model - fluids etc.
bgcMatl;

% Initialize parameter structure with constants, override as necessary.
prm = bgcConst();

% Initial conditions and other things related to the solver.
prm.solver.zo = 0; % [m]
prm.solver.zto = 0; % [m/s]
prm.solver.thetao = NaN;  % [K] not implemented.
prm.solver.tend = 3600*8; % s

% The background water column profile.
%prm.profile = bgcProfile(prm.const,'isopycnal/isothermal');
%prm.profile =  bgcProfile(prm.const,'pycnoclinic/isohaline');
%prm.profile =  bgcProfile(prm.const,'pycnoclinic/positive-upward haloclinic');
%prm.profile =  bgcProfile(prm.const,'pycnoclinic/positive-downward haloclinic');
%prm.profile = bgcProfile(prm.const,'levitus82World');
% ARGO profiles from Agulhas Current.
% these may need some pre-processing.  They may have areas of anomalously steep or weak gradients that
% confuse the output.
%argo = load('~/Dropbox/jakuba/veh/rafos_compressee/profiles/ArgoProfileMatrices_forMJ.mat');
%argo.ii = 56; % Profiles 56 and 83 decently bracket the variability.
%prm.profile = bgcProfile(prm.const,argo.argoPS(:,argo.ii),argo.argoTE(:,argo.ii),argo.argoPR(:,argo.ii));
% simplified Agulhas Current profile.
prm.profile = bgcProfile(prm.const,[35 35.5 35.25 34.5 34.75]',[25 17 13 8 3]',[0 200 400 800 2000]');

% The float as a whole has:
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
prm.h = 1; % [kg]
prm.CDs = 0; % [-] @@@@@
prm.CDf = 0.001; % [-]
prm.CDd = 0.1; % [-]  
prm.CDu = 0.1; % [-]
prm.D = 6*2.54/100; % [m]
prm.L = 1.5; % [m]
prm.As = pi*(prm.D/2)^2; % [m^2]
prm.Af = pi*prm.D*prm.L; % [m^2]


% Each component has:
% * a name 
% * mass, m [kg]
% * displacement, V [m^3], at 1 atm (zero if internal)
% * coefficient of thermal expansion, alpha [1/K]
% * an adiabatic compressibility, chi [1/Pa]
% * specific heat at constant pressur, cp [J/kg/K] (not used presently)
%

% reference values for water at the target pressure
[dens_kgpm3,T_K,P_Pa,~,~,~,S_PSU,~]=bgcInSitu(BALLAST_DEPTH,prm.profile);
chi = 1/sw_seck(S_PSU,T_K-273.15,P_Pa*1e-4)/1e5; % isopycnal at target depth

% Detailed component breakdown.
c = bgcInitComponent('Float');  % This is everything except the near-surface compressee.
c.V = (prm.D/2)^2*prm.L;
c.m = c.V*1000; % [kg] inclusive of all internal components.  Slightly positive.
c.rho = c.m/c.V;
c.chi = 0.1*1/aluminum.bulkModulus;  % @@@ arbitrary increased compressibility to account for elastic deformation of housing.
c.alpha = aluminum.coeffThermalExpansion;  % [m/K] Assumes internal components exert negligible internal pressure from their own expansion.  Certainly true in for gas-filled housing.
prm.components = bgcAddComponent(c);

% lockout volume (static volume of gas after lockout).
% compressibility of the structure surrounding this volume is assumed to
% have negligible effect on lockout volume.  @@@ valid?  Stable neutral depth is insensitive to this.
c = bgcInitComponent('Accumulator');
% @@@@@@@@@@ working here ^^^^
c.V = 0.1*(prm.D/2)^2*prm.L; % @@@@@@@@ this somehow has to be equal to the lockout volume.
c.m = c.V*1000; % [kg] inclusive of all internal components.  Slightly positive.
c.rho = c.m/c.V;
c.chi = 0.1*1/aluminum.bulkModulus;  % @@@ arbitrary increased compressibility to account for elastic deformation of housing.
c.alpha = aluminum.coeffThermalExpansion;  % [m/K] Assumes internal components exert negligible internal pressure from their own expansion.  Certainly true in for gas-filled housing.
prm.components = bgcAddComponent(c);


% Ballast. No effect on the above so long as ballast needs to be added and is
% added internal to the float.  This is assumed to apply beyond lockout depth.
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.

Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(T_K-prm.theta),(P_Pa-prm.const.atm)); % Volume of the system at depth.
Zc = prm.const.g*prmc.m - Vc*dens_kgpm3*prm.const.g; % (N) buoyancy (<0 indicates system is positive, >0 float is negative).
assert(Zc < 0,sprintf('Vehicle is negative at %.1f m (%.3f N).  This would require external volume and violate design assumptions.  Abort.',BALLAST_DEPTH,Zc));
fprintf(1,'Vehicle is positive.  Approx. %.1f kg margin ballast yields neutral at %.1f m\n',-Zc/prm.const.g,BALLAST_DEPTH);
f = bgcInitComponent('Internal ballast');
f.rho = inf;
f.alpha = 0;  % inside housing
f.chi = 0; % inside housing
f.m = -Zc/prm.const.g; % this will add mass only.
f.V = 0; % inside housing
prm.components = bgcAddComponent(f,prm.components);

% gas volumes not included in above.
c = bgcInitComponent('Air volume');  % 
c.V = 5000*CC2M3;  % [m^3] volume at 1 atm, T = prm.theta.
c.rho = 1.0; % [kg/m^3]  density at STP.
c.m = c.rho*c.V;
c.eos = @bgcVolumeIdealGas;
c.active = 1;
c.eventf = @bgcEventNone;
c.event_prm = {prm.profile,4000*CC2M3}; % {profile,lockoutVolume}
assert(c.V > c.event_prm{2}); 
prm.components = bgcAddComponent(c,prm.components);  % gas components handled separately.

%
% a discharge_rate property exists that could be used for semi-passive air-backed valve
%

% Descent controller for TR is just to thrust down until a certain depth is reached.
c = bgcInitComponent('descentController');
c.active = 0;  % 0 indicates active for this component only.  See bgcF.m
c.eventf = @bgcEventThrustDown;
c.event_prm = {CUTOFF_DEPTH,100}; % {stop depth [m],down thrust [N]}
prm.components = bgcAddComponent(c,prm.components);

% Drop weight not used for floats but needs to appear in component list and
% drop depth is used as first sample depth (bgcF.m)
c = bgcInitComponent('drop weight');
c.active = 0;
c.eventf = @bgcEventNone;
c.event_prm = {sampleDepths(1)};
prm.components = bgcAddComponent(c,prm.components);

% Surface and seafloor.  For convenience these are massless components
c = bgcInitComponent('bounds');
c.eventf = @bgcEventBounds;
c.event_prm = {0,1000.0};
prm.components = bgcAddComponent(c,prm.components);

% Controller must appear as last component.  For passive ballast does nothing except wait for
% sample timeout.
c = bgcInitComponent('controller');
c.active = 1;
c.eventf = @bgcEventNone; % always active
Zmax = 100; % [N] @@@ parameter passing issues btwn bgcF and feedback functions. 
c.event_prm = {sampleDepths,sampleDepthTol,sampleTime,sampleTimeLockout,@bgcFeedbackNone,NaN,NaN,NaN,NaN,Zmax,NaN,NaN};
prm.components = bgcAddComponent(c,prm.components);

prm
