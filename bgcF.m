function [yt,Zbuoyancy,Zdrag,Zthrust,rho,theta,p,m,Vf,thetaf,alpha,chi,cp,Re,zg,S] = bgcF(t,y,prm)
% Float dynamics, for use with ODE45.
%
% See notes, 2012-Dec-27 for derivations.
%
% Revision History
% 2012-12-27    mvj    Created.
% 2013-01-02    mvj    Modified to enable recovery of
%                      internal state after sim.
% 2022-03-14    mvj    add S output for isopycnal RAFOS study.
% 2025-03-31    mvj    add support for ideal gas
% 2025-05-07    mvj    add writing of log file.




% keep track of mission sample depths.
persistent tStartSample;
persistent izg;
persistent tlast;
if isempty(izg)
  izg = 1;
end

% Decompose state vector
zt = y(1);
z = y(2);
ize = y(3);
izte = y(4);


% Compute in situ properties.
% Displace the water column if specified.
if isfield(prm.profile,'displacement')
    switch prm.profile.displacement.type
      case 'sinusoidal'
        % not implemented because the dynamics do not account for motion of the ambient water column.
      case 'step'
        % get insitu S,T for displaced water.
        if t > prm.profile.displacement.step.tstep
            dz = prm.profile.displacement.step.step;
            [~,T_K,~,~,~,~,S,~] = bgcInSitu(z+dz,prm.profile);
            % move displaced water to current depth and compute insitu properties.
            p_db = sw_pres(z,0); 
            T_C = T_K - 273.15;
            T_C = T_C + dz*sw_adtg(S,T_C,p_db); % water will be adiabatically cooled or heated. 1st order approx.
            rho = sw_dens(S,T_C,p_db);
            p = p_db*1e4 + prm.const.atm;
            theta = T_C + 273.15;
        else
            [rho,theta,p,~,~,~,S] = bgcInSitu(z,prm.profile);
        end
    end
else
    [rho,theta,p,~,~,~,S] = bgcInSitu(z,prm.profile);
end

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


% Compute float temperature.
% @@@ No dynamics yet.  This assumes instantaneous equilibration.
thetaf = theta;

% Compute total volume and mass
[Vf,m] = deal(0); 
for c = 1:length(prm.components)

    cc = prm.components(c);
    fVolume = cc.eos;

    Vf = Vf + fVolume(cc.V,cc.alpha,cc.chi,(thetaf-prm.theta),(p-prm.const.atm),...
                      prm.theta,prm.const.atm,izg,cc.event_prm(:));    
    m = m + cc.m;
end

% 2025-03-31 bulk parameters for linear equation of state no long supported because
%            incompatible with ideal gas EoS.
[alpha,chi,cp] = deal(NaN);

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
bounds = prm.components(end-1);
assert(strcmp('bounds',bounds.name), ...
   'Bounds event does not appear as 2nd from last component in components list!');
descentCntrl = prm.components(end-3);
assert(strcmp('descentController',descentCntrl.name), ...
   'Descent controller does not appear as 3rd from last component in components list!');
cntrl = prm.components(end);
assert(strcmp('controller',cntrl.name), ...
    'Controller does not appear as last component in components list!');
dropweight = prm.components(end-2);
assert(strcmp('drop weight',dropweight.name), ...
    'Did not find drop weight at expected position in components list!');

% Default integral error rates.
ze = 0; zte = 0; 

% @@@ 2019/04/10 19:58:41  bgcIntegrator is deprecated.  There is no way presently to reset the integration
% now performed by ode45.  It can be saturated (see bgcFeedbackPIV).  Doesn't seem worth figuring out how to
% do this because the PIV scheme is clearly superior and avoids the need for any controller switching.
if ~descentCntrl.active
  % @@@ will probably have to alter this too.  This thing is set up here to engage and disengage the depth controller
  % @@@ once within a band, but we are writing a controller now that is always engaged.  might be able to handle that
  % @@@ using the active flag.
  ZthrustDescent = descentCntrl.event_prm{2};
  Zthrust = ZthrustDescent;
  zg = NaN;

elseif cntrl.active

  zFilter = cntrl.event_prm{1};
  zTol = cntrl.event_prm{2};
  ztTol = cntrl.event_prm{3};
  tSampleTime = cntrl.event_prm{4};
  Zmax = cntrl.event_prm{11};
  Zdead = cntrl.event_prm{14};

  if izg > length(zFilter)
    % hardcoded ascent.
    zg = NaN;
  else
    zg = zFilter(izg);
  end

  % Handle integrator windup.  Makes no sense while transiting.
  % @@@ not if PIV.
  %if abs(z-zg) > zTol
  %  bgcIntegrator(t,0,[],[]); % reset integrator.  This is critical.  Unclear if better than integral windup.
  %end

  if zg < 0 % interpret this as indicating closed loop ascent to surface.
    if bounds.active % @@@ does not differentiate between bottom and surface.
      izg = izg + 1;  % no guarantee this will be called only once.
      fprintf(1,'Arrived at surface.  Incrementing goal.\n');
    end
  end
      
  if isnan(zg)
      Zthrust = -Zmax;  % hard-coded ascent.
  else

      % bgcEventFilter starts and stops sampling when bounds on z, zt are met and z is close enough to
      % an element of zFilter.  It has no way of keeping track of which sample depth should be active
      % (izg is not part of the state).  bgcEventFilter can be used to shut the controller on and off
      % (other events too).  This block is still needed to time samples, though it does still violate
      % the scheme behind event-based interruptions.  Seems to be benign, at least relative to
      % toggling the controller here as opposed to using events.
    if isempty(tStartSample) && abs(z-zg) <= zTol && abs(zt) <= ztTol
      tStartSample = t;
      fprintf(1,'Starting new sample at depth: %.1f\n',zg);
    elseif (t-tStartSample) - tSampleTime  > 0 % sample done.
      tStartSample = [];
      izg = izg + 1;
      fprintf(1,'Sample complete.  New goal sample depth: %.1f\n',zFilter(izg));
    end
    
    % Run the controller.  
    fFeedback = cntrl.event_prm{6};
    [Zthrust,ze,zte] = fFeedback(t,zt,z,ize,izte,zg,cntrl.event_prm(7:end));
    % Simulate thruster deadband.  Controller could be written to compensate for this.
    if abs(Zthrust) < Zdead
        Zthrust = 0;
    end
    
  end
  
elseif dropweight.active % This has nothing to do with dropweight - it is the open loop thrust up between when the controller is engaged.
  %@@@@bgcIntegrator(t,0,[],[]); % reset integrator
  Zmax = cntrl.event_prm{10};
  Zthrust = -Zmax;
else
  %@@@@@bgcIntegrator(t,0,[],[]); % reset integrator
  Zthrust = 0;
  zg = NaN; % no goal depth if there is no controller running.
end

% Compute acceleration.  Zbouyancy is based on the profile and may include a free surface.
ztt = 1/(prm.h + m)*(Zbuoyancy + Zdrag + Zthrust);

% Create output vector including integral error terms.  
yt = [ztt; zt; ze; zte];

% Status.
fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%12.1f s %9.3f m %9.3f N  ',t,z,Zthrust);

% Write log.  Not insignificant slowdown.
fid = fopen('/tmp/cliosim.txt','a');
% [ztt; zt; ze; zte];
fprintf(fid,['%f ' ...
             '%f %f %f %f ' ...
             '%f %f %f ' ...
             '%f %f %f ' ...
             '%f %f %f ' ...
             '%f %f %f ' ...
             '%f %f %f ' ...
             '%f %f %f\n'], ...
        t,...
        zt,z,ize,izte, ...
        ztt,ze,zte, ...
        Zbuoyancy,Zdrag,Zthrust, ...
	rho,theta,p, ...
	m,Vf,thetaf, ...
	alpha,chi,cp, ...
	Re,zg,S);
fclose(fid);
