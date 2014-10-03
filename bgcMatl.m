% Some useful material properties.  Note the thermal expansion coefficients below are
% linear expansion coefficients.  We assume that the volumetric expansion coefficient
% is given by three times these values.
% Seawater properties derived from MATLAB seawater toolbox, which is now obsolete and should 
% be replaced with OSW Oceanographic Toolbox, http://www.teos-10.org/software.htm
% water (fresh): http://www.engineeringtoolbox.com/bulk-modulus-elasticity-d_585.html
% Aluminum: Wikipedia
% UHMW: Plastics International, TIVAR 1000.
% mild steel: http://www.ezlok.com/TechnicalInfo/MPCarbonSteel.html
% PFA: Dupont
% ESS syntactics (AZ grade): Compressive moduli provided - unclear if these are bulk modulus, K
% or uniaxial elastic modulus, E.  
% See TN600-5.doc - suggests most bouyancy foams have Poisson's Ratio \mu ~= 0.333, so that E = K.
% Syntactic AZ-34: ESS, 6000 m foam, AZ-Data-Sheet-Final-Draft-7-91.pdf
% Syntactic AZ-38: ESS, 7620 m foam, AZ_GRADE_TB.pdf  (2014/09/08 no longer available?)
% Syntactic HZ-34: ESS 7000 m foam, HZ-Data-Sheet-Final-7-10.pdf
% Syntactic Eccofloat DS-33, TrelleBorg, 6100 m foam.  Data sheet suggests material does not 
% follow the typical relationship for elastic behavior between the elastic modulus E, Poisson's 
% ratio \mu, and bulk modulus K: K = E/(3*(1-2*\mu)). \mu < 1/3 should result in E>K, but the 
% datasheet indicates \mu = 0.325 and E<K.
% See also http://www.synfoam.com/  for more foam options.
% water refers to fresh water, not seawater.  
water.coeffThermalExpansion = 69e-6; % [m/m/K]
water.bulkModulus = 2.15e9; % [Pa]  
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
% 2014/09/08  MVJ  Updated from new ESS datasheets, dated 7/11/2014.
syntacticAZ34.coeffThermalExpansion = 0; % @@@ unknown.
syntacticAZ34.bulkModulus = 2.75e9; % [Pa]  Assuming datasheet lists E and that \mu=.333
syntacticAZ34.density = 0.55e3; % [kg/m^3]
syntacticAZ34.thermalConductivity = NaN;
syntacticAZ34.specificHeat = NaN;
syntacticHZ34.coeffThermalExpansion = 0; % @@@ unknown.
syntacticHZ34.bulkModulus = 3.55e9; % [Pa]  Assuming datasheet lists E and that \mu=.333
syntacticHZ34.density = 0.55e3; % [kg/m^3]
syntacticHZ34.thermalConductivity = NaN;
syntacticHZ34.specificHeat = NaN;
% 2014/09/08  MVJ  AZ38 no longer available?
%syntacticAZ38.coeffThermalExpansion = 0; % @@@ unknown.
%syntacticAZ38.bulkModulus = 2.72e9; % [Pa]  Assuming datasheet lists E and that \mu=.333
%syntacticAZ38.density = 0.61e3; % [kg/m^3]
%syntacticAZ38.thermalConductivity = NaN;
%syntacticAZ38.specificHeat = NaN;
syntacticEccofloatDS33.coeffThermalExpansion = 0; % @@@ unknown.
syntacticEccofloatDS33.bulkModulus = 390000*PSI2PA; % [Pa] 2014/09/08  MVJ  Changed to K from E.
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
