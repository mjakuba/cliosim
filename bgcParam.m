function [prm,components] = bgcParam
% BGC Profiler parameters.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Added sampling parameters.  Needs revision to reflect real vehicle structure.


% Mission parameters
dropDepth = 6000;
sampleDepths = [20 120 320 520 720 920:500:5920]; % [m]
sampleDepthTol = 5; % +/- [m]
sampleDepthRateTol = 0.1; % +/- [m/s]
sampleTime = 2250; % [s]
sampleTimeLockout = 50; % [s]

% Conversions
PSI2PA = 6894.75729;
LB2KG = 0.453592;
FT2M = 12*2.54/100;
IN2M = 2.54/100;
DEGF2K =  5/9;

% Constants
prm.const.g = 9.81; % [m/s^2]
prm.const.nu = 1.83e-6; % [m^2/s, 0 C]
prm.const.mu = 1.88e-3; % [Pa s, 0 C]
prm.const.atm = 101325; % [Pa] atmospheric pressure

% Some useful material properties.  Note the thermal expansion coefficients below are
% linear expansion coefficients.  We assume that the volumetric expansion coefficient
% is given by three times these values.
% Aluminum: Wikipedia
% UHMW: Plastics International, TIVAR 1000.
% mild steel: http://www.ezlok.com/TechnicalInfo/MPCarbonSteel.html
% PFA: Dupont
% Syntactic AZ-34: Esyntactic, 6230 m foam.  Compressive modulus provided, but not bulk modulus.
% Syntactic AZ-38: Esyntactic, 7620 m foam.  Compressive modulus provided, but not bulk modulus.
% Syntactic Eccofloat DS-33, TrelleBorg, 6100 m foam.
water.coeffThermalExpansion = 69e-6; % [m/m/K]
water.bulkModulus = 2.2e9; % [Pa] 
water.density = 999; % [kg/m^3]
water.thermalConductivity = 0.58; % [W/m/K]
water.specificHeat = 4187; % [J/kg/K]
rnalater.density = 1340; % [kg/m^3]
aluminum.coeffThermalExpansion = 23.1e-6; % [m/m/K @ 25 C]
aluminum.bulkModulus = 76e9; % [Pa]
aluminum.density = 2.7e-3*100^3; % [kg/m^3 @ 20 C]
aluminum.thermalConductivity = 237; % [W/m/K]
mildSteel.coeffThermalExpansion = 11.5e-6; % [m/m/K @ 25 C]
mildSteel.bulkModulus = 140e9; % [Pa]
mildSteel.density = 7.870e-3*100^3; % [kg/m^3 @ 20 C]
mildSteel.thermalConductivity = 51.9; % [W/m/K]
uhmw.coeffThermalExpansion = 0.00011/DEGF2K; % [m/m/K @ 23 C]
uhmw.bulkModulus = 77750*PSI2PA; 
uhmw.density = 58.01*LB2KG/FT2M^3; % [kg/m^3 @ 23 C]
uhmw.thermalConductivity = NaN;
syntacticAZ34.coeffThermalExpansion = 0; % @@@ unknown.
syntacticAZ34.bulkModulus = 2.62e9; % [Pa]
syntacticAZ34.density = 0.55e3; % [kg/m^3]
syntacticAZ34.thermalConductivity = NaN;
syntacticAZ34.specificHeat = NaN;
syntacticAZ38.coeffThermalExpansion = 0; % @@@ unknown.
syntacticAZ38.bulkModulus = NaN; % [Pa]
syntacticAZ38.density = 0.61e3; % [kg/m^3]
syntacticAZ38.thermalConductivity = NaN;
syntacticAZ38.specificHeat = NaN;
syntacticEccofloatDS33.coeffThermalExpansion = 0; % @@@ unknown.
syntacticEccofloatDS33.bulkModulus = 375000*PSI2PA; % [Pa]
syntacticEccofloatDS33.density = 0.50e3; % [kg/m^3]
syntacticEccofloatDS33.thermalConductivity = NaN;
syntacticEccofloatDS33.specificHeat = NaN;
pfa.coeffThermalExpansion = 13e-5; % [m/m/K]
pfa.bulkModulus = 480e6/3/(1-2*0.46); % [Pa] Poisson's ratio for PTFE used.
pfa.density = 2150; % [kg/m^3]
pfa.thermalConductivity = 0.195; % [W/m/K]
pfa.specificHeat = 1172; % [J/kg/K]
polycarbonate.coeffThermalExpansion = 70e-6; % [m/m/K]
polycarbonate.bulkModulus = 3.1e9; % [Pa]
polycarbonate.density = 1200; % [kg/m^3]
polycarbonate.thermalConductivity = 0.2; % [W/m/K]
polycarbonate.specificHeat = 1200; % [J/kg/K]
pvc.coeffThermalExpansion = 50.4e-6; % [m/m/K]
pvc.bulkModulus = 2.41e9/3/(1-2*0.3825); % [Pa] 
pvc.density = 1300; % [kg/m^3]
pvc.thermalConductivity = 0.147; % [W/m/K]
pvc.specificHeat = 1355; % [J/kg/K]
ti2.coeffThermalExpansion = 8.6e-6; % [m/m/K]
ti2.bulkModulus = 1.05e11/3/(1-2*0.33); % [Pa] 
ti2.density = 4510; % [kg/m^3]
ti2.thermalConductivity = 21.79; % [W/m/K]
ti2.specificHeat = 500; % [J/kg/K]
hdpe.coeffThermalExpansion = 200e-6; % [m/m/K]
hdpe.bulkModulus = 1.07e09/3/(1-2*0.4101); % [Pa] 
hdpe.density = 952; % [kg/m^3]
hdpe.thermalConductivity = 0.461; % [W/m/K]
hdpe.specificHeat = 1796; % [J/kg/K]
mineraloil.coeffThermalExpansion = 1/3*6.4e-4; % [m/m/K]
mineraloil.bulkModulus = 1.9e9; % [Pa] 
mineraloil.density = 875; % [kg/m^3]
mineraloil.thermalConductivity = 0.1; % [W/m/K]
mineraloil.specificHeat = 0.45*water.specificHeat; % [J/kg/K]

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
prm.CDf = 0.01; % [-]
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
c = bgcInitComponent('SUPR 1');
c.m = 4; 
c.V = 0.0023;
c.rho = c.m/c.V;
c.alpha = polycarbonate.coeffThermalExpansion;
c.chi = 1/polycarbonate.bulkModulus;
c.cp = polycarbonate.specificHeat;
prm.components = bgcAddComponent(c);
c.name = 'SUPR 2';
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('Filter Stack 1');
c.m = 25; 
c.V = 0.019;
c.rho = c.m/c.V;
c.alpha = polycarbonate.coeffThermalExpansion;
c.chi = 1/polycarbonate.bulkModulus;
c.cp = polycarbonate.specificHeat;
prm.components = bgcAddComponent(c);
c.name = 'Filter Stack 2';
prm.components = bgcAddComponent(c,prm.components);

% Battery packs are internal, so no volume as far as simulation 
% is concerned.  This is for Oceanserver long-format battery packs.
% Each is 0.095 kWh.  6 packs are necessary for filtering (500 Wh).
% Assume for now that 500 Wh for propulsion is reasonable.
c = bgcInitComponent('Battery Pack');
c.m = 12*1.425*LB2KG;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('Battery Housing 1');
c.m = 10; 
c.V = (31*IN2M)*pi*(3.5/2*IN2M)^2;
c.rho = c.m/c.V;
c.alpha = ti2.coeffThermalExpansion;
c.chi = 1/ti2.bulkModulus;  % @@@ actual housing will be considerably less stiff than a solid rod of Ti.
c.cp = NaN;
prm.components = bgcAddComponent(c,prm.components);
c.name = 'Battery Housing 2';
prm.components = bgcAddComponent(c,prm.components);
c.name = 'Battery Housing 3';
prm.components = bgcAddComponent(c,prm.components);
c.name = 'Battery Housing 4';
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('Shell');
c.m = 20; 
c.V = 0.0084;
c.rho = c.m/c.V;
c.alpha = pfa.coeffThermalExpansion;
c.chi = 1/pfa.bulkModulus;
c.cp = pfa.specificHeat;
prm.components = bgcAddComponent(c,prm.components);

% @@@ structure seems a significant underestimate.
% @@@ Perhaps a fair amount of Ti used here?
% @@@ solid model dubious.  Made this Ti.  
c = bgcInitComponent('Structure (Placeholder)');
c.V = 0.01;
c.rho = ti2.density;
c.m = c.V*c.rho;
c.alpha = ti2.coeffThermalExpansion;
c.chi = 1/ti2.bulkModulus;
c.cp = ti2.specificHeat;
prm.components = bgcAddComponent(c,prm.components);

% @@@ Several other smaller housings not included.

% Add RNA Later bag.
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
f = bgcInitComponent('flotation');
f.rho = syntacticAZ34.density;
f.alpha = syntacticAZ34.coeffThermalExpansion;
f.chi = 1/syntacticAZ34.bulkModulus;
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);
f.V = 0.186; % [m^3]
% bgcNeutral(6000.0,prm.profile,prmc,f)*1.02;  % 2% reserve buoyancy at depth.
f.m = f.rho*f.V;
%114.380126; % [kg]
prm.components = bgcAddComponent(f,prm.components);

% Thrust for in the case of no drop weight.
c = bgcInitComponent('descentController');
c.eventf = @bgcEventThrustDown;
ZthrustDescent = 50; % [N] @@@ needs to be checked - this is in forward direction with non-zero advance velocity.
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
