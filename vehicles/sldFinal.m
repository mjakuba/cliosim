function [prm,components] = bgcParam
% BGC Profiler parameters.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Added sampling parameters.  Needs revision to reflect real vehicle structure.
% 2014-10-03    mvj    Using this to suck in material components and associated properties from Solidworks.
% 2014-10-08    mvj    Hard-coded solidworks output file here instead of in bgcAddSldComponents.m
% 2015-08-31    mvj    Updated drag coefficient based on 1/4-scale model testing.  And updated to something
%                      close to what appeared in the AGU 2014 Fall meeting concept.

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
% * skin friction coefficient, CDf [-], referenced to surface area
% * downcast quadratic drag coefficient, CDd [-], referenced to frontal area
% * upcast quadratic drag coefficient, CDu [-], referenced to frontal area
% * surface area, As [m^2]
% * frontal area, Af [m^2]
% * characteristic diameter, D [m]
prm.theta = 300; % [K]
prm.h = 10; % [kg]
prm.CDf = 0.001; % [-]
prm.As = 0.0; % [m^2] Surface area.  Set to zero because model test drag coefficients include frictional losses.
%prm.Af = 0.22; % [m^2]  Base area.
prm.Af = 0.83; % [m^2] Frontal area.  Seems that based on model testing we're getting a base drag that suggests
               % separation as soon as the body begins to taper.  There is no value to the taper on the afterbody.
	       % (There is of course considerable value on the forebody).
prm.D = 0.5; % [m] @@@ sets Stokes drag - presumably based on sphere or something?  Needs to be updated in any case,
             % but should not matter much for energy calcs.
prm.CDd = 0.3; % [-]  Drag coefficient referenced to frontal area. (based on 1/4-scale model tests).  This includes
               % any frictional component.
prm.CDu = 0.3; % [-]  As for CDd.


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

% H-beam concept more or less as at Mechanical Design Review, 2015-12-16.
% 2015/12/15 21:21:07  @@@ what other components are missing - chassies in particular?  Also need to 
% mirror SUPRs and various mounts.
!cp /home/jakuba/Dropbox/Clio/vehicle/sld/tmp/H-beam_concept20160321163145.sldtxt /tmp/tmp.sldtxt # Trelleborg EL-34
%!cp /home/jakuba/Dropbox/Clio/vehicle/sld/tmp/H-beam_concept20160321164025.sldtxt /tmp/tmp.sldtxt # Trelleborg TG-32
!sed -i s#Delrin\ 2700\ NC010,#Delrin\ 2700\ NC010\ # /tmp/tmp.sldtxt
!sed -i s#Transducer\ Backing\ Plate,\ ITC-3013#Transducer\ Backing\ Plate\ ITC-3013# /tmp/tmp.sldtxt
prm.components = ...
    bgcAddSldComponents('/tmp/tmp.sldtxt', ...
    prm.components); 


% 2015/09/02 17:53:53  Battery packs not in solid model.  Adding them here for now.  
% 2016/03/21 20:43:15  Update for battery modules.
c = bgcInitComponent('Battery Pack');
c.m = 4*2*4.33+2*0.44+2*1.78;  % Current concept will be able to fit up to 32*2 packs.
prm.components = bgcAddComponent(c,prm.components);

% entrained water?  Probably significant - we are likely 
% to haul a significant amount of surface water around 
% with us.  However, simulation scheme does not support 
% entraining water of a specific density.  To do so it would 
% need to not use the events scheme (easy) and maintain some state
% between simulation segments (doable)


% Add status print.
c = bgcInitComponent('Status Printer');
c.eventf = @bgcEventPrintStatus;
c.event_prm = {};
prm.components = bgcAddComponent(c,prm.components);


% Add RNA Later bag.
c = bgcInitComponent('RNAlater');
c.rho = rnalater.density;
%c.V = 10/1000; % 10 l based on other SUPR configurations.
c.V = 0/1000; % 2015/09/02 13:25:29  Current design calls for trivial quantity of RNALater.
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
c.V = 0/1000*(rnalater.density-1030)/(1030-c.rho); % 2015/09/02 13:25:29  Current design calls for trivial quantity of RNALater. 
c.m = c.V*c.rho;
c.alpha = water.coeffThermalExpansion;
c.chi = 1/water.bulkModulus;
c.cp = water.specificHeat;
c.discharge_rate = c.V/(sampleTime*length(sampleDepths));
c.eventf = @bgcEventFilter;
c.event_prm = {sampleDepths,sampleDepthTol,sampleDepthRateTol,sampleTime,sampleTimeLockout};
prm.components = bgcAddComponent(c,prm.components);


% 2014-10-22  Feedback at the conceptual stage we really need is just how positive or negative the current design is
%             Much easier to add margin ballast than margin flotation.  Note though that RNA later discharge and
%             compensator are still in here as from proposal - that's 110 kg of mass, but with minimal effect on
%             net buoyancy.
BALLAST_DEPTH = 100;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
[rho,theta,p] = bgcInSitu(BALLAST_DEPTH,prm.profile); % seawater at 3000 m - for concept studies, we'd like the float to be
                                             % neutral at half depth.
Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); % Volume of the float at depth.
Zc = prm.const.g*prmc.m - Vc*rho*prm.const.g; % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
if Zc < 0
  fprintf(1,'Vehicle is positive.  Approx. %.1f kg margin ballast yields neutral at %.1f m\n',-Zc/prm.const.g,BALLAST_DEPTH);
else
  fprintf(1,'Vehicle is negative at %.1f m.  %.1f N additional floatation required.\n',BALLAST_DEPTH,Zc);
end

% 2015/09/03 21:41:41 Added this note so I dont wonder later why this doesn't work!
warning('this has insufficient thrust to overcome buoyancy at depth and negative at surface.  Sim takes forever.');

% 2014/10/22 20:39:14 Don't need the rest of this for conceptual design for neutral.
fprintf(1,'Paused.  Hit any key to add/rm floatation for neutral at %.1f m and continue.\n',BALLAST_DEPTH);
pause;
%return

% Add ballast/flotation as necessary to make vehicle neutral at 3000 m.
SYNTACTIC = syntacticTrelleborgEL34;
if Zc < 0
  f = bgcInitComponent('ballast (negatively buoyant)');
  f.rho = aluminum.density;
  f.alpha = aluminum.coeffThermalExpansion;
  f.chi = 1/aluminum.bulkModulus;
  margin = 0.995;
else
  f = bgcInitComponent('ballast (positively buoyant)');
  f.rho = SYNTACTIC.density;
  f.alpha = SYNTACTIC.coeffThermalExpansion;
  f.chi = 1/SYNTACTIC.bulkModulus;
  margin = 1.005;
end
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
f.V = bgcNeutral(BALLAST_DEPTH,prm.profile,prmc,f)*margin;
f.m = f.rho*f.V;
prm.components = bgcAddComponent(f,prm.components);

% 2016/03/22 14:49:43  MVJ  Determine Si Oil volume required to make Clio roughly isopycnal.
% 2016/03/22 19:04:11  MVJ  This is still pretty far off from isopycnal because of thermal expansion
%                           and because compressibility of seawater varies with temperature
%                           and pressure.  Presumably there is some optimal depth we should use the
%                           seawater bulk modulus for.  Note that this will yield a positively buoyanct
%                           vehicle because we're recomputing flotation for the slight positive buoyancy of
%                           silicone oil.
% 2016/03/23 20:36:00  MVJ  Silicone oil is highly responsive to temperature.  A postive vehicle at the surface
%                           may be very negative at 20 m.   Have to use a fairly shallow ballast depth.
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.
fprintf('Vehicle''s net bulk modulus and net thermal coefficient of expansion, w/o compressee: kappa = %.2g; alpha = %.2g\n',1/prmc.chi,prmc.alpha);
% 2018/10/17 12:26:30. Data from clio005 says the actual compressibility was about prmc.chi = 1/4.8e9;  This requires
% huge amounts of silicone oil (340 kg) to compensate.  Makes sense, should be about half vehicle mass.  Prediction
% here was a bulk modulus of 2.5e9, much closer to seawater.  Unclear why we were so far off.   The data from
% Trelleborg actually indicates a slightly lower bulk modulus than this simulation assumed.  This model has only 2 SUPRs.
o = bgcInitComponent('silicone oil compressee');
o.rho = siliconeoil.density;
o.alpha = siliconeoil.coeffThermalExpansion;
o.chi = 1/siliconeoil.bulkModulus;
o.V = prmc.V*(1/seawater.bulkModulus - prmc.chi)/(1/siliconeoil.bulkModulus - 1/seawater.bulkModulus);
o.m = o.V*siliconeoil.density;
fprintf(1,'Silicone oil added to make vehicle isopycnal: %.1f cc, %.1f kg.  Hit any key to continue.\n', ...
    o.V*100^3,o.m);
pause
prm.components = bgcAddComponent(o,prm.components);


% Thrust for in the case of no drop weight.
c = bgcInitComponent('descentController');
c.eventf = @bgcEventThrustDown;
% Overriding this now to get vehicle down.  Auto-flotation computation may make it neutral at 6000 m, but it
% is wicked buoyant at the surface.  (This is a totally bogus vehicle so no issue yet.)
% 2015/08/31 19:42:08  Reset to 50 N.  Vehicle design is bogus, but ballasting for neutral at 3000 m is probably a
% generally better idea than ballasting for margin at 6000 m from an energy perspective since we're not likely to
% get anywhere close to isopycnal.  (Emergency Drop weight will have to be sized to compensate).
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
%Zmax = 50; % [N] might be generous.
Zmax = 500; % definitely generous - just to get sim to run.
c.event_prm = {sampleDepths,sampleDepthTol,sampleTime+50,sampleTimeLockout,Kp,Kd,Ki,Zff,Zmax};
prm.components = bgcAddComponent(c,prm.components);

