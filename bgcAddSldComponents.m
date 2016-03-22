function components = bgcAddSldComponents(fmacro,components)
% bgcAddSldComponents(fmacro)
%
% fmacro 
%
% 2014/10/03 18:18:39  Read these in from a the output of the sldwrks macro.
% 2014/10/06 20:43:17  Screwed up archived copy of first run with bogus vehicle.  Need to restore from dbox.  Eh.
% who cares?
% 2014-10-07    mvj    Functionalized.  
% 2014-10-08    mvj    Computation or extraction of custom bulk modulus now dealt with in solidworks macro.


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
    % 2014-10-08    mvj    Seems like there are precision problems on some parts.  I do think the part is being
    %                      updated temporarily if necessary.  This showed up on a part that did not require updating.
    %                      Changed assert to a warning.
    %assert(abs(sld.mass(n)/sld.volume(n) - sld.density(n)) < eps(1e5)); 
    if abs(sld.mass(n)/sld.volume(n) - sld.density(n)) > 0.1
      warning(sprintf('%s: mass/volume differs from material density by %.16f kg/m^3', ...
	  sld.componentName{n}, ...
	  sld.mass(n)/sld.volume(n) - sld.density(n)))
    end
    

    c.m = sld.mass(n);
    c.V = sld.volume(n);
    c.rho = sld.density(n);
    c.alpha = sld.coeffThermalExpansion(n);
   
    % Bulk modulus.  Since there are various schemes in play for using effective/approximate bulk moduli
    % at least check that bulk modulus is not 0.  That definitely indicates a missing material or effective property
    % and the simulation will fail.  Alternately could set compressibility to 0 in this case and issue a warning.
    if (sld.bulkModulus(n) == 0)
      %error(sprintf('Bulk Modulus for solid component %s is zero.',c.name{:}));
      warning(sprintf('Bulk Modulus for solid component %s is zero!  Setting compressibility to zero.',c.name{:}));
      c.chi = 0.0;
    else
      c.chi = 1/sld.bulkModulus(n);
    end
  
    c.cp = sld.specificHeat(n);
    
  elseif strcmp(sld.componentType(n),'internal')
    
    c.m = sld.mass(n);
    % All other values ignored and left at defaults.

  elseif strcmp(sld.componentType(n),'OEM')  % external OEM part/assy.
    
    c.m = sld.mass(n);
    c.V = sld.volume(n);
    % All other values ignored and left at defaults.

  elseif strcmp(sld.componentType(n),'housing')  % Compressibility 
    
    c.m = sld.mass(n);
    c.V = sld.volume(n);
    c.chi = 1/sld.bulkModulus(n);
    % All other values ignored and left at defaults.
    
  else
   error('Unsupported component type.');
 end
   
 components = bgcAddComponent(c,components);

end
