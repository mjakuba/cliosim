function [prm] = bgcParam
% Armada SBIR phase 1 float concept.
%
% Focus is on mesoscale.
%
% Revision History
% 2025-03-31    mvj    Created.
% 2025-04-06    mvj    This version uses PIV and executes a complete cycle back to the surface.

% o ascent to 0 isn't right from energy use perspective.  Want to ascend to above unstable neutral depth,
%   but not to settle at that depth; rather to drive through.  Driving all the way to the surface is probably close enough.
% x tuning is way too slow, but seems to work well.  Was able to use piv_design to choose gains.
% o need to turn off controller in some way, probably use a deadband on thrust, which is realistic anyway.  alternatively
%   a deadband on depth and depth rate errors.

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
prm.solver.tend = 3300; % s  Stop sim before final ascent or it slows way down.


% Mission parameters.  Adapt Clio mission to float mission.  Relevant logic is in bgcF.m
% Ascent is triggered after all samples are complete (bgcF.m).
BALLAST_DEPTH = 500; % I think this should be deeper than any target depth?  Yes, this will be the completely compressed depth.
NEUTRAL_DEPTH = 50;
CUTOFF_DEPTH = 5; % [m] initial descent to below unstable near-surface neutral depth
sampleDepths = [400 -10 250 -10]; % [m]  0s make this run very slow.  unclear why, maybe in bgcEventFilter?
sampleDepthTol = 0.1; % +/- [m] Sample timer starts once within this band.
sampleDepthRateTol = 0.1; % +/- [m/s] and depth rate below this figure.  If too small, replay won't catch this.
sampleTime = 1000; %7200; % [s] time to remain at sample depth.
sampleTimeLockout = 60; % [s] ?
COMPRESSEE_VOLUME = 100*CC2M3; % [m^3] Maximum volume change.


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
%prm.profile = bgcProfile(prm.const,[35 35.5 35.25 34.5 34.75]',[25 17 13 8 3]',[0 200 400 800 2000]');
prm.profile = bgcProfile(prm.const,'isopycnal/isothermal');  % this will induce only pressure effects on device.  Seawater in situ density will be affected primarily by pressure, with small effect from salinity.

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
prm.theta = 300; % [K]
prm.h = 1; % [kg]
prm.CDs = 0; % [-] 
prm.CDf = 0.001; % [-] 
prm.CDd = 0.1; % [-] 
prm.CDu = 0.1; % [-]
prm.D = 6*2.54/100; % [m]
prm.L = 1.5; % [m]
prm.As = pi*prm.D*prm.L; % [m^2]
prm.Af = pi*(prm.D/2)^2; % [m^2]

% Each component has:
% * a name 
% * mass, m [kg]
% * displacement, V [m^3], at 1 atm (zero if internal)
% * coefficient of thermal expansion, alpha [1/K]
% * an adiabatic compressibility, chi [1/Pa]
% * specific heat at constant pressur, cp [J/kg/K] (not used presently)
%

% reference values for water at the target pressure
[dens_kgpm3,T_K,P_Pa,~,~,~,S_PSU,~]=bgcInSitu(NEUTRAL_DEPTH,prm.profile);
chi = 1/sw_seck(S_PSU,T_K-273.15,P_Pa*1e-4)/1e5; % isopycnal at target depth

% Detailed component breakdown.
c = bgcInitComponent('Float');  % This is everything except the near-surface compressee.
c.V = (prm.D/2)^2*prm.L;
c.m = c.V*1000; % [kg] inclusive of all internal components.  Slightly positive.
c.rho = c.m/c.V;
% If chi is very large, device is unstable about BALLAST_DEPTH and returns to the surface.  
c.chi = 20*1/aluminum.bulkModulus;  % @@@ factor of 10 increase is what was computed for Clio housings.  Shallow will be softer.
c.alpha = aluminum.coeffThermalExpansion;  % [m/K] 
c.eos = @bgcVolumeLinear;
prm.components = bgcAddComponent(c);

% Ballast for unstable equilibrium at neutral depth.
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
Vc = bgcVolumeLinear(prmc.V,prmc.alpha,prmc.chi,(T_K-prm.theta),(P_Pa-prm.const.atm)); % Volume of the system at neutral depth.
Vg = bgcVolumeIdealGas(COMPRESSEE_VOLUME,NaN,NaN,(T_K-prm.theta),(P_Pa-prm.const.atm),prm.theta,prm.const.atm); % Volume of the compressee at neutral depth.
Zc = prm.const.g*prmc.m - (Vc+Vg)*dens_kgpm3*prm.const.g; % (N) buoyancy (<0 indicates system is positive, >0 float is negative).
assert(Zc < 0,sprintf('Vehicle is negative at %.1f m (%.3f N).  This would require external volume and violate design assumptions.  Abort.',BALLAST_DEPTH,Zc));
fprintf(1,'Vehicle is positive.  Approx. %.1f kg margin ballast yields neutral at %.1f m\n',-Zc/prm.const.g,BALLAST_DEPTH);
f = bgcInitComponent('Internal ballast');
f.rho = inf;
f.alpha = 0;  % inside housing
f.chi = 0; % inside housing
f.m = -Zc/prm.const.g; % this will add mass only.
f.V = 0; % inside housing
f.eos = @bgcVolumeLinear;
prm.components = bgcAddComponent(f,prm.components);

% Now compute required lockout volumes for each of the target depths.
clear lockoutVolumes;
for n=1:length(sampleDepths)
    
    % seawater properties at target depth
    [dens_kgpm3,T_K,P_Pa,~,~,~,S_PSU,~]=bgcInSitu(sampleDepths(n),prm.profile);

    % system at depth exclusive of compressee
    [prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
    Vc = bgcVolumeLinear(prmc.V,prmc.alpha,prmc.chi,(T_K-prm.theta),(P_Pa-prm.const.atm)); % Volume of the system at depth.
    Zc = prm.const.g*prmc.m - Vc*dens_kgpm3*prm.const.g; % (N) buoyancy (<0 indicates system is positive, >0 float is negative).
    assert(Zc > 0,sprintf('Vehicle is positive at %.1f m (%.3f N).  Cannot add compressee volume.  Abort.',sampleDepths(n),Zc));

    % assuming air is massless
    lockoutVolumes(n) = Zc/(dens_kgpm3*prm.const.g);

    % make it such that lockout volume changes after each surfacing.
    % this overwrites depth changes without surfacing.  No deep volume changes.
    if sampleDepths(n) <= 0 
        lockoutVolumes(n) = lockoutVolumes(n-1);
    end

end

%dump
%return


% add flotation in the form of residual compressee volume.
% Modeled as massless.
c = bgcInitComponent(sprintf('Compressee'));  
c.V = COMPRESSEE_VOLUME;  % [m^3] volume at 1 atm, T = prm.theta.
c.rho = 1; % [kg/m^3]  density at STP.
c.m = 0; % set to zero such that ballasting above is accurate. 
c.eos = @bgcVolumeIdealGas;  % The valve is effectively in here - use the lockout volume parameter
c.eventf = @bgcEventNone; % @@@ something will need to advance the goal lockout volume.  that should happen in bgcF.m
c.event_prm = {prm.profile,lockoutVolumes}; %{prm.profile,lockoutVolumes};
assert(all(c.V > c.event_prm{2})); 
prm.components = bgcAddComponent(c,prm.components); 


%
% a discharge_rate property exists that could be used for semi-passive air-backed valve
%

% Descent controller for TR is just to thrust down until a certain depth is reached.
% @@@ possibly might want to reengage this to deal with large initial air volume?
c = bgcInitComponent('descentController');
c.active = 1;  % 1 means disabled for tihs component only.  Use PIV instead.
prm.components = bgcAddComponent(c,prm.components);

% Drop weight not used for floats but needs to appear in component list.
c = bgcInitComponent('drop weight');
c.active = 0;
c.eventf = @bgcEventNone;
c.event_prm = {NaN};  % NaN skips using the dropweight depth for initial descent.
prm.components = bgcAddComponent(c,prm.components);

% Surface and seafloor.  For convenience these are massless components
c = bgcInitComponent('bounds');
c.eventf = @bgcEventBounds;
c.event_prm = {0,1000.0,false,true}; % {surf, bot, stop-sim-on-surfacing, stop-sim-on-grounding}
prm.components = bgcAddComponent(c,prm.components);

% Controller must appear as last component.  For passive ballast does nothing except wait for
% sample timeout.
c = bgcInitComponent('controller');
c.active = 1;
%c.eventf = @bgcEventNone; % always active
c.eventf = @bgcEventFilter;  % this doesn't work.
Zmax = 15; % [N] 
Zballast = 0; % [N]  This is the assumed and is very nearly correct below cutoff.
Kp = 0.1; Kv = 10; Ki = 0.7;  % parameters from piv_design_armada.m
ztmax = 1.0; % [m/s]
Zww = -0.5*1000*(prm.CDf*prm.As + prm.CDu*prm.Af); % Not used in bgcFeedbackPIV.m  Apparently destabilizing.
Zdead = 0.3; %0.1; %1.0; % [N] hard to imagine producing less than 0.1 N thrust consistently.  Probably more like 1 N.
c.event_prm = {sampleDepths,sampleDepthTol,sampleDepthRateTol,sampleTime,sampleTimeLockout,@bgcFeedbackPIV,Kp,Kv,Ki,ztmax,Zmax,Zballast,Zww,Zdead};
prm.components = bgcAddComponent(c,prm.components);

dump

