function [prm,components] = bgcParam
% isobaric Rossby float, idealized
%
% Revision History
% 2022-03-02    mvj    Created.

% o use profile from Heather for off Madagascar, but first use an isothermal and isohaline profile
% o use some kind of controller to get 

% Mission parameters
BALLAST_DEPTH = 400; % HF 2022-06-06: 260-300 dbar or 560-640 dbar

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

% reference values for water at the target pressure
[dens_kgpm3,T_K,P_Pa,~,~,~,S_PSU,~]=bgcInSitu(BALLAST_DEPTH,prm.profile);
chi = 1/sw_seck(S_PSU,T_K-273.15,P_Pa*1e-4)/1e5; % isopycnal at target depth

% compressee design.
% 2022-07-27 Get some parameters from the root workspace to allow iteration.
ENSEMBLE = true;
if ENSEMBLE
    CF = evalin('base','CF');
else
    CF = 2.956e-6*1e-4; % [1/Pa] from CICESE ballast sheets (and used in Seascan design, 2022-06-01)
end
VF = 12125*1e-6; % [m^3] 2022-06-01 design
% This is what you would do, not what the design necessarily is.
%MARGIN = 1.0;
%CS = MARGIN*chi;  % design compressibility
%VC = 990*1e-6;  % 2022-06-01 design; This is mostly the structure. the active volume is small.
%comp = CS*(VC+VF)-CF*VF;  % for a piston this is 1/K*(pi*D^2/4)^2  [m^5/N] or [m^3/Pa] or [cc/Pa]*1e-6
%CC = comp/VC;

% @@@@@@@@@@@@@@@@@ Wrong!  We need CF from Seascan design + variability, but CS should be
% and this is then still assuming that the compressee components have no variability, which is probably
% a pretty bad assumption given the use of plastics
% also need to update to 7075.

% The float, excluding the compressee and final ballasting to make it neutral.
% @@@ masses and volumes a little arbitrary from CICESE ballast sheets
%c = bgcInitComponent('Float');
%c.V = VF;
%c.m = 11217*1e-3;
%c.rho = c.m/c.V;
%c.chi = CF;
%c.alpha = 1/3*0.096e-4; % this is for glass tube, from Rossby et al., 1985.  Waiting on HF for new data.
%prm.components = bgcAddComponent(c);

% compressee - lumped approximate version.
%c = bgcInitComponent('Compressee');
%c.V = VC; % 2022-07-27  This works out if VC an CC are computed as above, but this is not the active volume alone.
%c.rho = 0.5*aluminum.density; % guess at adjustment for oil volume. 
%c.m = c.rho*c.V;
%c.chi = CC; % this assumes the active volume change completely dominates any contraction of the rest of the device.
%c.alpha = aluminum.coeffThermalExpansion; % may actually be dominated by oil?  oil ~ 10 more thermally responsive.
%prm.components = bgcAddComponent(c,prm.components);

% Detailed component breakdown.
c = bgcInitComponent('DropWeightHolder');  % this is the spacer, was UHMW, now Al
c.V = 133*1e-6;
c.rho = aluminum.density;
c.m = c.rho*c.V;
c.chi = 1/aluminum.bulkModulus;
c.alpha = aluminum.coeffThermalExpansion;
prm.components = bgcAddComponent(c);

c = bgcInitComponent('EndCap');
c.V = 192*1e-6;
c.rho = aluminum.density;
c.m = c.rho*c.V;
c.chi = 1/aluminum.bulkModulus;
c.alpha = aluminum.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('Float');  % This is the glass and everything inside but does not include the endcap.
cn = {'DropWeightHolder','EndCap'};
c.V = VF; % [m^3]
for n=1:length(cn)
    ii = length(prm.components);
    assert(~isempty(strfind(cn,prm.components(ii).name)));
    c.V = c.V - prm.components(ii).V;
end
c.m = 12.6-2.6 % [kg] 2022-06-01 design, total estimated weight - estimated compressee weight.
c.rho = c.m/c.V;
c.chi = CF; % defined above.
c.alpha = borosilicate_glass.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('CompresseeActiveVolume') % components are handled below.  This captures the action of the spring.
SPRINGK = 554.4e2; % N/m  2022-06-01 design.
DPISTON = 2.051e-2; % m  2022-06-01 design.
LPISTON = 4.2e-2; % m
c.V = LPISTON*pi*DPISTON^2/4; % total throw. 
c.rho = 0;
c.m = 0; % oil is modeled separately.
c.chi = (1/SPRINGK*(pi*DPISTON^2/4)^2)/c.V; 
c.alpha = 0;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('CompresseeDiaphram');
c.V = 2.3*1e-6;
c.rho = rubber.density;
c.m = c.V*c.rho;
c.chi = 1/rubber.bulkModulus;
c.alpha = rubber.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('CompresseeOilTube');
c.V = 7*1e-6;
c.rho = polyurethane.density;
c.m = c.rho*c.V;
c.chi = 1/polyurethane.bulkModulus;
c.alpha = polyurethane.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('CompresseeOilTubeFitting');
c.V = 2.36*1e-6;
c.rho = polypropylene.density;
c.m = c.rho*c.V;
c.chi = 1/polypropylene.bulkModulus;
c.alpha = polypropylene.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);


c = bgcInitComponent('CompresseeOilTubePlug');
c.V = 1*1e-6;
c.rho = polypropylene.density;
c.m = c.rho*c.V;
c.chi = 1/polypropylene.bulkModulus;
c.alpha = polypropylene.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);


% compressee oil.  From 2114 Compressee Detailed Design, 2022-03-16, when device is bottomed out,
% the assembly contains 3.72 cc + 4.5 cc of oil.  At max extension, it contains 11.12 cc + 4.5 cc of oil.
% The compensation tube is 1.4" ID tube, ~ 7" long.  When full its volume would be 5.6 cc.  That's actually not
% enough since the delta between full extension and bottoming out is 7.4 cc
% Take at face value for now.
% 2022-06-01 update.  Volume of oil is 3.4 cc (oil in the rod) + 13.6 (oil below the bellofram when fully relaxed) + 9 (minimum compensator oil) + an unk margin volume in the comp.  Perhaps 30 cc total?
c = bgcInitComponent('CompresseeOil');
c.V = 30*1e-6;
c.rho = mineraloil.density;
c.m = c.rho*c.V;
c.chi = 1/mineraloil.bulkModulus;
c.alpha = mineraloil.coeffThermalExpansion;  % linear coeff (see bgcMatl.m)
prm.components = bgcAddComponent(c,prm.components);

c = bgcInitComponent('CompresseeAl6061Assy'); % all the Al parts, mostly 6061 plus the spring 7075.
cn = {'CompresseeOil','CompresseeDiaphram','CompresseeOilTubePlug','CompresseeOilTubeFitting','CompresseeOilTube'};
c.V = 990e-6; % 2022-06-01 design
for n=1:length(cn)
    ii = length(prm.components);
    assert(~isempty(strfind(cn,prm.components(ii).name)));
    c.V = c.V - prm.components(ii).V;
end
c.rho = aluminum.density;
c.m = c.rho*c.V;
c.chi = 1/aluminum.bulkModulus;
c.alpha = aluminum.coeffThermalExpansion;
prm.components = bgcAddComponent(c,prm.components);

%c = bgcInitComponent('CompresseeAir');
% @@@@@@@@@ do we need to do this?

% Ballast. No effect on the above so long as ballast needs to be added and is added internal to the float.
% @@@@ assume ballast is installed internal to the floats.  This may not be true.
prmc = prm;
[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.

% 2022/04/01 12:24:06 results - all expressed as volumetric.
% HHF got net alpha = 1.3537e-5 /C (volumetric)
% w oil, I get 1.7706e-5
% w/o oil, I get 1.67e-5

%dump; return;

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
