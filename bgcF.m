function [yt,Zbuoyancy,Zdrag,Zthrust,rho,theta,p,m,Vf,thetaf,alpha,chi,cp,Re,zg] = bgcF(t,y,prm)
% Float dynamics, for use with ODE45.
%
% See notes, 2012-Dec-27 for derivations.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Modified to enable recovery of
%                      internal state after sim.

% keep track of mission sample depths.
persistent tStartSample;
persistent izg;
persistent tlast;
if isempty(izg)
  izg = 0;  % marks descent.
  tlast = 0;
end

% Decompose state vector
zt = y(1);
z = y(2);

% Compute in situ properties.
[rho,theta,p] = bgcInSitu(z,prm.profile);

% Modify the mass and volume of any active components.  This supports
% discharging only.
% @@@ discharge rate applies to surface volume, not volume at depth, 
% @@@ but that is almost certainly insignificant.
for c = 1:length(prm.components)
  if prm.components(c).active && prm.components(c).discharge_rate ~= 0
    dt = (t-prm.components(c).activate_time);  % time since last activation.
    
    prm.components(c).V = prm.components(c).V - dt*prm.components(c).discharge_rate;
    if prm.components(c).V < 0;
      prm.components(c).V = 0;
    end
    prm.components(c).m = prm.components(c).rho*prm.components(c).V;
    
  end
end

% Derive bulk quantities for profiler.
[m,V,alpha,chi,cp] = ...
    bgcBulkParam(prm.components);

% Compute float temperature.
% @@@ No dynamics yet.  This assumes instantaneous equilibration.
thetaf = theta;

% Compute volume at this temperature and pressure
Vf = bgcVolume(V,alpha,chi,(thetaf-prm.theta),(p-prm.const.atm));

% Buoyancy force.
Zbuoyancy = m*prm.const.g - Vf*rho*prm.const.g;

% Drag force.
Re = prm.D*abs(zt)/prm.const.nu;
ZdragStokes = -6*pi*prm.const.mu*prm.D/2*zt;
q = 0.5*rho*zt*abs(zt);
if zt > 0 % downcast
  ZdragTurb = -(q*prm.CDf*prm.As + q*prm.CDd*prm.Af);
else % upcast
  ZdragTurb = -(q*prm.CDf*prm.As + q*prm.CDu*prm.Af);
end
if Re <= 1
  Zdrag = ZdragStokes;
elseif Re < 1e3 % @@@ bogus model for transitional Re.
  Zdrag = ZdragStokes + ZdragTurb;
else
  Zdrag = ZdragTurb;
end

% Control.
% @@@ reference to specific components is really ugly.  Consider redoing components struct
descentCntrl = prm.components(end-3);
assert(strcmp('descentController',descentCntrl.name), ...
   'Descent controller does not appear as 3rd from last component in components list!');
cntrl = prm.components(end);
assert(strcmp('controller',cntrl.name), ...
    'Controller does not appear as last component in components list!');
dropweight = prm.components(end-2);
assert(strcmp('drop weight',dropweight.name), ...
    'Did not find drop weight at expected position in components list!');

if ~descentCntrl.active
  % @@@ will probably have to alter this too.  This thing is set up here to engage and disengage the depth controller
  % @@@ once within a band, but we are writing a controller now that is always engaged.  might be able to handle that
  % @@@ using the active flag.
  ZthrustDescent = descentCntrl.event_prm{2};
  Zthrust = ZthrustDescent;
elseif cntrl.active

  zFilter = cntrl.event_prm{1};
  zTol = cntrl.event_prm{2};
  tSampleTime = cntrl.event_prm{3};
  Zmax = cntrl.event_prm{10};
  
  if izg == 0
    zg = dropweight.event_prm{1}; % initial descent depth.
  elseif izg > length(zFilter)
    % hardcoded ascent.
    %zg = -1000;  weird results somehow sim loops back on itself.
    zg = NaN;
  else
    zg = zFilter(izg);
  end

  % Handle integrator windup.  Makes no sense while transiting.
  % @@@ not if PIV.
  %if abs(z-zg) > zTol
  %  bgcIntegrator(t,0,[],[]); % reset integrator.  This is critical.  Unclear if better than integral windup.
  %end
  
  
  if isnan(zg)
    Zthrust = -Zmax;
  else

    % Start sample timer once within depth band.
    if isempty(tStartSample) && abs(z-zg) <= zTol
      tStartSample = t;
    elseif (t-tStartSample) - tSampleTime  > 0 % sample done.
      tStartSample = [];
      izg = izg + 1;
      fprintf(1,'New goal sample depth: %.1f\n',zg);
    end
    
    fFeedback = cntrl.event_prm{5};
    Zthrust = fFeedback(t,zt,z,zg,cntrl.event_prm(6:end));
    
  end
  
elseif dropweight.active % This has nothing to do with dropweight - it is the open loop thrust up between when the controller is engaged.
  bgcIntegrator(t,0,[],[]); % reset integrator
  Zmax = cntrl.event_prm{10};
  Zthrust = -Zmax;
else
  bgcIntegrator(t,0,[],[]); % reset integrator
  Zthrust = 0;
end

% Compute acceleration.
ztt = 1/(prm.h + m)*(Zbuoyancy + Zdrag + Zthrust);

% Create output vector.  
yt = [ztt; zt];
