function bgcPlotBuoyancy(prm)
warning('Unfinished!')
% really though to look at the contributions you need to separate the density changes in the water from just pressure
% and temperature effects.

Z = [0:6000];
THETA = [0:35]+271;

[prmc.m,prmc.V,prmc.alpha,prmc.chi,prmc.cp] = bgcBulkParam(prm.components);  % compute effective parameters.

figure(1); clf reset;

% first hold theta constant at the minimum
theta = 270;
for n = 1:length(Z)
  z = Z(n);
  
  [rho,~,p] = bgcInSitu(z,prm.profile);
  
  % Volume of the float at depth.
  Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); 
   % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
  Zc(n) = prm.const.g*prmc.m - Vc*rho*prm.const.g;
end

plot(Zc,Z); hold on;

% now hold theta constant at the maximum.
theta = 310;
for n = 1:length(Z)
  z = Z(n);
  
  [rho,~,p] = bgcInSitu(z,prm.profile);
  
  % Volume of the float at depth.
  Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); 
   % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
  Zc(n) = prm.const.g*prmc.m - Vc*rho*prm.const.g;
end

plot(Zc,Z); hold on;


% both effects (use Levitus profile)
for n = 1:length(Z)

  z = Z(n);
  [rho,theta,p] = bgcInSitu(z,prm.profile);
  
  % Volume of the float at depth.
  Vc = bgcVolume(prmc.V,prmc.alpha,prmc.chi,(theta-prm.theta),(p-prm.const.atm)); 
   % (N) buoyancy (<0 indicates float is positive, >0 float is negative). 
  Zc(n) = prm.const.g*prmc.m - Vc*rho*prm.const.g;
end

plot(Zc,Z);


legend('270 K','310 K','Levitus Profile');
grid on;
xlabel('Z (N), Negative is Buoyant');
ylabel('Depth (m)');
