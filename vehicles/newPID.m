function [prm,components] = bgcParam
% BGC Profiler parameters.
%
% Revision History
% 2019-04-09    mvj    Create to explore inner/outer loop depth control.



% Mission parameters
dropDepth = 2010; % used by OL descent.
sampleDepths = [2000 1500 1000 600 400 315 200 100 30]; % [m]
sampleDepthTol = 5; % +/- [m]
sampleDepthRateTol = 0.1; % +/- [m/s]
sampleTime = 1920; % [s]
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
prm.CDd = 0.3; % [-]
prm.CDu = 0.3; % [-]
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

% hard-coded filename within for the moment.
prm.components = bgcInitComponent('Null Component to initialize prm.components');  % @@@ bogus...

c = bgcInitComponent('Lumped vehicle');
c.rho = 1050;
c.V = 633/1000;
c.m = c.V*c.rho;
%c.alpha = 0; % no info
c.chi = 1/3.6e9; % approximate, from data.
c.alpha = 8.9e-5;
%c.chi = 8.1e9; % from data but seems very high.  2017-11-01  No idea what this number is.  Should be O(1e-10).
c.cp = NaN;
prm.components = bgcAddComponent(c,prm.components);

% compressee
c = bgcInitComponent('Silicone Oil');
c.rho = siliconeoil.density;
c.V = 200/1000;  % roughly 200 l gives a net bulk modulus of 2.2 GPa (water) for k = 3.6e9.
c.m = c.V*c.rho;
c.alpha = siliconeoil.coeffThermalExpansion; % no info
c.chi = 1/siliconeoil.bulkModulus; % crude - this varies a lot.
c.cp = siliconeoil.specificHeat;
prm.components = bgcAddComponent(c,prm.components);


% ballasting.  This should be small relative to lumped vehicle above. 
% Add ballast/flotation as necessary to make vehicle neutral at 3000 m.
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


% Active descent.
c = bgcInitComponent('descentController');
c.eventf = @bgcEventThrustDown;
ZthrustDescent = 300;  % N
c.event_prm =  {dropDepth,ZthrustDescent};
prm.components = bgcAddComponent(c,prm.components);

% Drop weight is not included in computation of flotation.
% Not used - m=0.
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
c.event_prm = {0,2300.0};
prm.components = bgcAddComponent(c,prm.components);

% Controller is always last component.  This is a component mostly
% for convenience - it has no mass or volume.
c = bgcInitComponent('controller');
c.eventf = @bgcEventControl;
Kp = 100; Kd = 50; Ki = 0.01; 
Zff = 150; % [N]
Zmax = 500; % [N]
c.event_prm = {sampleDepths,sampleDepthTol,sampleTime+50,sampleTimeLockout,@bgcFeedbackPID,Kp,Kd,Ki,Zff,Zmax};
prm.components = bgcAddComponent(c,prm.components);

% issue is bgcF parses out Kp,Ki and then calls bgcFeedback to do PID.  There is presently no way to call a
% generic controller from within.  instead it always calls bgcFeedback (PID)
% also bgcEventControl turns the controller on/off when within a depth band instead of incrementing the goal 
% after a timeout.

