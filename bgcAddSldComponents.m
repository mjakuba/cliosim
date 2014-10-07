% 2014/10/03 18:18:39  Read these in from a the output of the sldwrks macro.
% 2014/10/06 20:43:17  Screwed up archived copy of first run with bogus vehicle.  Need to restore from dbox.  Eh.
% who cares?
%fmacro = '/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/CLIO-199-0000.SLDASM.txt';  % hard-coded for now
%fmacro = '/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/proposalVehicle.txt';
fmacro = '/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/proposal_vehicle20140710124610.sldtxt';
%autoread('/home/jakuba/Dropbox/Clio/vehicle/sld/tmp/CLIO-199-0000.SLDASM.txt');
% autoread won't work right now because of scientific notation and also doubles being interpreted as %d if the first
% line happens to have spat out an integer.
[sld.componentName,sld.componentType,sld.mass,sld.volume,sld.matlName,sld.density,sld.elasticModulus,sld.poissonsRatio, ...
      sld.shearModulus,sld.bulkModulus,sld.coeffThermalExpansion,sld.thermalConductivity,sld.specificHeat] = textread(fmacro, ...
    '%[^,],%[^,],%f,%f,%[^,],%f,%f,%f,%f,%f,%f,%f,%f','headerlines',2);

% Add components from solidworks.  Check for some errors.
% This supports only bulk materials 2014/10/03 20:51:51
% Will need to add assertions for internal parts (will show up in type field).
for n=1:length(sld.componentName)
  
  c = bgcInitComponent(sld.componentName(n));  
  if strcmp(sld.componentType(n),'solid')

    % I think I have this fixed, but in some cases seems like bodies (as opposed to parts) had custom material properties
    % assigned and this catches those.
    assert(abs(sld.mass(n)/sld.volume(n) - sld.density(n)) < eps(1e5)); 

    c.m = sld.mass(n);
    c.V = sld.volume(n);
    c.rho = sld.density(n);
    c.alpha = sld.coeffThermalExpansion(n);
   
    % Deal with bulk modulus.
    if ~isnan(sld.bulkModulus(n))
      warning('Using custom bulk modulus for solid component.');
      c.chi = 1/c.bulkModulus(n);  % @@@ not sure this is a good idea - at least warn that we are using a hard-coded
      % custom property for some component as delivered.
    else
      bulkModulus = sld.elasticModulus(n)/(3*(1-2*sld.poissonsRatio(n)));
      c.chi = 1/bulkModulus;
    end
  
     c.cp = sld.specificHeat(n);
    
  elseif strcmp(sld.componentType(n),'internal')
    
    c.m = sld.mass(n);
    % All other values ignored and left at defaults.
   
  else
   error('Unsupported component type.');
 end
   
  
   
   
   % Cheating here - functionalize this later.
   prm.components = bgcAddComponent(c,prm.components);
   
  
end
