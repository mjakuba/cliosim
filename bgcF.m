function [yt,Zbuoyancy,Zdrag,Zthrust,rho,theta,p,m,Vf,thetaf,alpha,chi,cp,Re] = bgcF(t,y,prm)
% Float dynamics, for use with ODE45.
%
% See notes, 2012-Dec-27 for derivations.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Modified to enable recovery of
%                      internal state after sim.

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
  if prm.components(c).active
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
elseif Re > 1e3 % @@@ bogus model for transitional Re.
  Zdrag = ZdragStokes + ZdragTurb;
else
  Zdrag = ZdragTurb;
end

% Control.
% @@@ reference to specific components is really ugly.  Consider redoing components struct
% @@@ into struct with named fields.
cntrl = prm.components(end);
assert(strcmp('controller',cntrl.name), ...
    'Controller does not appear as last component in components list!');
dropweight = prm.components(end-2);
assert(strcmp('drop weight',dropweight.name), ...
    'Did not find drop weight at expected position in components list!');
if cntrl.active
  Kp = cntrl.event_prm{5};
  Kd = cntrl.event_prm{6};
  Ki = cntrl.event_prm{7};
  Zff = cntrl.event_prm{8};
  % Determine goal.  Assume it is closest to current state.
  % This relies on the controller being good enough to avoid
  % falling into the "well" of another goal depth.
  [nul,izg] = min(abs(cntrl.event_prm{1}-z)); 
  zg = cntrl.event_prm{1}(izg);
  Zthrust = bgcFeedback(t,zt,z,zg,Kp,Kd,Ki,Zff);
elseif dropweight.active
  bgcIntegrator(t,0,[],[]); % reset integrator
  Zmax = cntrl.event_prm{9};
  Zthrust = -Zmax;
else
  bgcIntegrator(t,0,[],[]); % reset integrator
  Zthrust = 0;
end

% Compute acceleration.
ztt = 1/(prm.h + m)*(Zbuoyancy + Zdrag + Zthrust);

% Create output vector.
yt = [ztt; zt];
