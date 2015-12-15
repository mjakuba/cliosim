% BGC Profiler simulation master file.
%
%
% Revision History
% 2012-12-27    mvj    Created.
% 2014-10-07    mvj    Fixed output event list wrt wierd extra non-terminal events sometimes generated 
%                      at start of new simulation segment


% Set up paths
addpath('./seawater');
addpath('./util');

% Load parameters
prm = bgcParam;

% Set initial conditions.
zo = 0; % [m]
zto = 0; % [m/s]
thetao = prm.theta; % [K]  Not used.

% Solver parameters
yo = [zto; zo];
tstart = 0; % [s]
tend = 50000; % [s]
% @@@ Event functions sometimes problematic for large maxstep.  Not sure how to improve that behavior without
% @@@ two event functions devoted to each depth interval of interest so that zero crossings are guaranteed.
odeOptions = odeset('MaxStep',2, ...
    'Events',@(t,y) bgcEvents(t,y,prm), ...
    'Stats','on');

% Solve.  Every event is terminal, solver is used repeatedly and output aggregated.
tout = tstart;
yout = yo';
teout = [];
yeout = [];
ieout = [];
stop = false;
[Zbuoyancy,Zdrag,Zthrust, ...
      rho,theta,p, ...
      mf,Vf,thetaf, ...
      alphaf,chif,cpf, ...
      Re] = deal(NaN);
tize_last = [];
ize_last = [];
clear functions;  % Persistent variables are used for some timing operations.
while tout(end) < tend && ~stop
  
  [t,y,te,ye,ie] = ode45(@(t,y) bgcF(t,y,prm),[tstart tend],yout(end,:),odeOptions);
  
  % Recover some non-state data from simulation for this segment.  Requires
  % careful treatment of integrator in feedback.
  [tize,ize] = bgcIntegrator(t(end),0); % get integral error at end of this segment
  bgcIntegrator(t(1),0,tize_last,ize_last); % reset integral error back to beginning of this segment
  ii = length(tout)+[1:length(t)];
  [Zbuoyancy(ii,1),Zdrag(ii,1),Zthrust(ii,1), ...
	rho(ii,1),theta(ii,1),p(ii,1), ...
	mf(ii,1),Vf(ii,1),thetaf(ii,1), ...
	alphaf(ii,1),chif(ii,1),cpf(ii,1), ...
	Re(ii,1)] = deal(NaN*ones(length(t),1));
  for n=1:length(t)
    nn = length(tout) + n;
    [nul,Zbuoyancy(nn),Zdrag(nn),Zthrust(nn), ...
	  rho(nn),theta(nn),p(nn), ...
	  mf(nn),Vf(nn),thetaf(nn), ...
	  alphaf(nn),chif(nn),cpf(nn), ...
	  Re(nn)] = ...
	bgcF(t(n),y(n,:),prm);
  end
  bgcIntegrator(t(end),0,tize,ize); % Necessary?
  tize_last = tize;
  ize_last = ize;
  
  % (De)Activate any components whose events occurred.  All events 
  % are assumed to be toggles.
  for iie = 1:length(ie)

    c = ie(iie);    
    % Some kind of wierd bug.  Often an event that just happened will
    % happen on the start of the next cycle, but ode45 won't terminate
    % on it but then reports it!  Ignore these and make sure they don't
    % show up in the output used for plotting.
    if t(end) - te(iie) > 1.0
      continue
    else
       teout = [teout; te(iie)];
       yeout = [yeout; ye(iie,:)];
       ieout = [ieout; ie(iie)];
     end
    
    % Permanently update state of active components.
    if prm.components(c).active
      dt = (te(iie)-prm.components(c).activate_time); 
      prm.components(c).V = prm.components(c).V - dt*prm.components(c).discharge_rate;
      if prm.components(c).V < 0;
	prm.components(c).V = 0;
      end
      prm.components(c).m = prm.components(c).rho*prm.components(c).V;
    end
    
    % End simulation upon impacting seafloor or reaching surface.
    if strcmp(prm.components(c).name,'bounds')
      stop = true;
    end
    
    % Toggle components that triggered event.
    prm.components(c).active = ~prm.components(c).active;
    prm.components(c).activate_time = te(iie);

    % Display event.
    if prm.components(c).active
      fprintf(1,'Event (t=%.1f, z=%.1f zt = %0.2f): %s Engaged (component %d).\n', ...
	  te(iie),ye(iie,2),ye(iie,1),prm.components(c).name,ie(iie));
    else
      fprintf(1,'Event (t=%.1f, z=%.1f zt = %0.2f): %s Disengaged (component %d).\n', ...
	  te(iie),ye(iie,2),ye(iie,1),prm.components(c).name,ie(iie));
    end
    
  end
  
  % Store simulation output.
  tout = [tout; t];
  yout = [yout; y];
  % teout = [teout; te];
  % yeout = [yeout; ye];
  % ieout = [ieout; ie];
  
  % Update initial time for next segment of sim.
  tstart = tout(end);
  
end

% Dump variables to workspace
dump;

% Plot results.
bgcPlot;
