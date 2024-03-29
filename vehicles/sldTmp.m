function [prm,components] = bgcParam
% BGC Profiler parameters.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Added sampling parameters.  Needs revision to reflect real vehicle structure.
% 2014-10-03    mvj    Using this to suck in material components and associated properties from Solidworks.
% 2014-10-08    mvj    Hard-coded solidworks output file here instead of in bgcAddSldComponents.m


% Mission parameters
dropDepth = 6000;
sampleDepths = [20 120 320 520 720 920:500:5920]; % [m]
sampleDepthTol = 5; % +/- [m]
sampleDepthRateTol = 0.1; % +/- [m/s]
sampleTime = 2250; % [s]
sampleTimeLockout = 50; % [s]

% Load conversion constants.
% @@@ shouldn't be necessary to use because solidworks spits out data in SI units by default.
bgcConversions;

% Still need this file for things that aren't in the model - fluids etc.
bgcMatl;

% Initialize parameter structure with constants, override as necessary.
prm = bgcConst();

% The background water column profile.
prm.profile = bgcProfile(prm.const);

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
prm.theta = 300; % [K]
prm.h = 10; % [kg]
prm.CDs = 0; % [-]
prm.CDf = 0.001; % [-]
prm.CDd = 0.1; % [-]
prm.CDu = 0.1; % [-]
prm.As = 1.0; % [m^2]
prm.Af = 0.22; % [m^2]
prm.D = 0.5; % [m]

% Each component has:
% * a name 
% * mass, m [kg]
% * displacement, V [m^3], at 1 atm (zero if internal)
% * coefficient of thermal expansion, alpha [1/K]
% * an adiabatic compressibility, chi [1/Pa]
% * specific heat at constant pressur, cp [J/kg/K] (not used presently)
%

% hard-coded filename within for the moment.
prm.components = bgcInitComponent('Null Component to initialize prm.components');  % @@@ bogus...
% 2014/10/22 18:40:29  Proposal vehicle uses derived flotation, which has been commented out.
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/proposal_vehicle20140710124610.sldtxt',prm.components); % proposal vehicle
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/20142210145919.sldtxt',prm.components); % skinny fatboy (120 N negative)
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/20142710092830.sldtxt',prm.components); % fat fatboy (800 N positive)
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/20142210134222.sldtxt',prm.components); % 4 pod.
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/placeholder_supr_assy20142210163259.sldtxt',prm.components); % one SUPR
%prm.components = bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/20142210164450.sldtxt',prm.components); % one filter holder
prm.components = ...
    bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/H-beam_concept20141031154827.sldtxt', ...
   prm.components); % H-beam concept without SUPRs, and too small mid-section.
% prm.components = ...
%     bgcAddSldComponents('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/placeholder_supr_assy20143110151700.sldtxt', ...
%     prm.components); % 9-filter SUPR, some missing components and motor still internal.

% % Battery packs are internal, so no volume as far as simulation 
% % is concerned.  This is for Oceanserver long-format battery packs.
% % Each is 0.095 kWh.  6 packs are necessary for filtering (500 Wh).
% % Assume for now that 500 Wh for propulsion is reasonable.
% c = bgcInitComponent('Battery Pack');
% c.m = 12*1.425*LB2KG;
% prm.components = bgcAddComponent(c,prm.components);

% c = bgcInitComponent('Battery Housing 1');
% c.m = 10; 
% c.V = (31*IN2M)*pi*(3.5/2*IN2M)^2;
% c.rho = c.m/c.V;
% c.alpha = ti2.coeffThermalExpansion;
% c.chi = 1/ti2.bulkModulus;  % @@@ actual housing will be considerably less stiff than a solid rod of Ti.
% c.cp = NaN;
% prm.components = bgcAddComponent(c,prm.components);
% c.name = 'Battery Housing 2';
% prm.components = bgcAddComponent(c,prm.components);
% c.name = 'Battery Housing 3';
% prm.components = bgcAddComponent(c,prm.components);
% c.name = 'Battery Housing 4';
% prm.components = bgcAddComponent(c,prm.components);



% Add RNA Later bag.
c = bgcInitComponent('RNAlater');
c.rho = rnalater.density;
%c.V = 10/1000; % 10 l based on other SUPR configurations.
c.V = 0/1000; % 10 l based on other SUPR configurations.
c.m = c.V*c.rho;
c.alpha = water.coeffThermalExpansion; % @@@ no info
c.chi = 1/water.bulkModulus; % @@@ no info
c.cp = NaN;
c.discharge_rate = c.V/(sampleTime*length(sampleDepths));
c.eventf = @bgcEventFilter;
c.event_prm = {sampleDepths,sampleDepthTol,sampleDepthRateTol,sampleTime,sampleTimeLockout};
prm.components = bgcAddComponent(c,prm.components);

% Add compensator for RNA Later discharge.  Otherwise vehicle gets too buoyant as it ascends.
% This could also just be DI water, which we will have to discharge anyway for water samples.
% @@@ DI water discharge is important - it will make the vehicle less buoyant and so is dangerous.
c = bgcInitComponent('Compensator');
c.rho = water.density;
%c.V = 10/1000*(rnalater.density-1030)/(1030-c.rho); 
c.V = 0/1000*(rnalater.density-1030)/(1030-c.rho); 
c.m = c.V*c.rho;
c.alpha = water.coeffThermalExpansion;
c.chi = 1/water.bulkModulus;
c.cp = water.specificHeat;
c.discharge_rate = c.V/(sampleTime*length(sampleDepths));
c.eventf = @bgcEventFilter;
c.event_prm = {sampleDepths,sampleDepthTol,sampleDepthRateTol,sampleTime,sampleTimeLockout};
prm.components = bgcAddComponent(c,prm.components);


% entrained water?  Probably significant - we are likely 
% to haul a significant amount of surface water around 
% with us.  However, simulation scheme does not support 
% entraining water of a specific density.  To do so it would 
% need to not use the events scheme (easy) and maintain some state
% between simulation segments (doable)

% Handle flotation as a derived quantity to yield a
% neutral float at desired depth.
% @@@ Leave this in here for now - it may have some utility as a feedback mechanism to solidworks design.
% f = bgcInitComponent('flotation');
% f.rho = syntacticEccofloatDS33.density;
% f.alpha = syntacticEccofloatDS33.coeffThermalExpansion;
% f.chi = 1/syntacticEccofloatDS33.bulkModulus;
% prmc = prm;
% [prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
% f.V = bgcNeutral(6000.0,prm.profile,prmc,f)*1.02;  % 2% reserve buoyancy at depth.
% f.m = f.rho*f.V;
% prm.components = bgcAddComponent(f,prm.components);
%
% 2014-10-22  Feedback at the conceptual stage we really need is just how positive or negative the current design is
%             Much easier to add margin ballast than margin flotation.  Note though that RNA later discharge and
%             compensator are still in here as from proposal - that's 110 kg of mass, but with minimal effect on
%             net buoyancy.
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
[rho,theta,p] = bgcInSitu(3000,prm.profile); % seawater at 3000 m - for concept studies, we'd like the float to be
                                             % neutral at half depth.
Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); % Volume of the float at depth.
Zc = prm.const.g*prmc.m - Vc*rho*prm.const.g; % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
if Zc < 0
  fprintf(1,'Vehicle is positive.  Approx. %.1f kg margin ballast yields neutral at 3000 m',-Zc/prm.const.g);
else
  warning(sprintf('Vehicle is negative at 3000 m.  %.1f N additional floatation required.',Zc));
end

keyboard

% 2014/10/22 20:39:14 Don't need the rest of this for conceptual design for neutral.
return

% Thrust for in the case of no drop weight.
c = bgcInitComponent('descentController');
c.eventf = @bgcEventThrustDown;
% Overriding this now to get vehicle down.  Auto-flotation computation may make it neutral at 6000 m, but it
% is wicked buoyant at the surface.  (This is a totally bogus vehicle so no issue yet.)
%ZthrustDescent = 50; % [N] @@@ needs to be checked - this is in forward direction with non-zero advance velocity.
ZthrustDescent = 500;
c.event_prm =  {dropDepth,ZthrustDescent};
prm.components = bgcAddComponent(c,prm.components);

% Drop weight is not included in computation of flotation.
c = bgcInitComponent('drop weight');
c.m = 0; % [kg]
c.rho = mildSteel.density;
c.V = c.m/c.rho;
c.alpha = mildSteel.coeffThermalExpansion;
c.chi = 1/mildSteel.bulkModulus;
c.cp = NaN;
c.discharge_rate = 10; % [m^3/s] effectively instantaneous
c.eventf = @bgcEventDropWeight;
c.event_prm = {dropDepth};
prm.components = bgcAddComponent(c,prm.components);

% Surface and seafloor.  For convenience these are massless components
c = bgcInitComponent('bounds');
c.eventf = @bgcEventBounds;
c.event_prm = {0,6500.0};
prm.components = bgcAddComponent(c,prm.components);

% Controller is always last component.  This is a component mostly
% for convenience - it has no mass or volume.
c = bgcInitComponent('controller');
c.eventf = @bgcEventControl;
Kp = 10; Kd = 50; Ki = 0.01; 
Zff = 0; % [N]
Zmax = 50; % [N] might be generous.
c.event_prm = {sampleDepths,sampleDepthTol,sampleTime+50,sampleTimeLockout,Kp,Kd,Ki,Zff,Zmax};
prm.components = bgcAddComponent(c,prm.components);

