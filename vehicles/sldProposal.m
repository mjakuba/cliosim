function [prm,components] = bgcParam
% BGC Profiler parameters.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Added sampling parameters.  Needs revision to reflect real vehicle structure.
% 2014-10-03    mvj    Using this to suck in material components and associated properties from Solidworks.
% 2014-10-06    mvj    Tried to model proposal vehicle in solidworks to see if can reproduce proposal results.
%                      Ran using foam computed by simulation, not foam in vehicle model.
% 2014-10-07    mvj    Oops - accidentally ran the thrusting down vehicle variant of the proposal vehicle.  Fixed
%                      to original proposal vehicle.


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
bgcAddSldComponents;

% Add RNA Later bag.  
% 2014/10/07 17:23:41  Noticed this is not modeled correctly - RNALater is discharged
%                      throughout filtering instead of just at the end.  Not particularly
%                      relevant to current design, and leaving this here anyway because
%                      this vehicle is supposed to show proposal results can be recreated
%                      using solidworks in the loop.
c = bgcInitComponent('RNAlater');
c.rho = rnalater.density;
c.V = 10/1000; % 10 l based on other SUPR configurations.
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
c.V = 10/1000*(rnalater.density-1030)/(1030-c.rho); 
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
f = bgcInitComponent('flotation');
f.rho = syntacticEccofloatDS33.density;
f.alpha = syntacticEccofloatDS33.coeffThermalExpansion;
f.chi = 1/syntacticEccofloatDS33.bulkModulus;
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
f.V = bgcNeutral(6000.0,prm.profile,prmc,f)*1.02;  % 2% reserve buoyancy at depth.
f.m = f.rho*f.V;
prm.components = bgcAddComponent(f,prm.components);

% Was not a component of proposed vehicle - it was to use a drop weight, but current
% version of code requires that a descent controller exists as a component immediately
% above the drop weight.  For proposal vehicle, just make associated thrust 0.
c = bgcInitComponent('descentController');
c.eventf = @bgcEventThrustDown;
ZthrustDescent = 0;
c.event_prm =  {dropDepth,ZthrustDescent};
prm.components = bgcAddComponent(c,prm.components);

% Drop weight is not included in computation of flotation.
% 2014/10/07 16:53:04  Altered back to proposal vehicle (m=20 instead of m=0).
c = bgcInitComponent('drop weight');
c.m = 20; % [kg]
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

