% BGC Profiler simulation master file.
%
%
% Revision History
% 2013-02-05    mvj    Created.


figure(1); clf reset;
subplot(211)
plot(tout,yout(:,2));
set(gca,'ydir','reverse');
tticklabel('abs',3600);
ylabel('Depth (m)');
grid on;

subplot(212)
plot(tout,yout(:,1));
set(gca,'ydir','reverse');
ylabel('Depth Rate (m/s)');
tticklabel('abs',3600);
grid on;

figure(2); clf reset;
subplot(211)
plot(tout,[Zbuoyancy,Zdrag,Zthrust]);
legend('Zbuoyancy','Zdrag','Zthrust');
tticklabel('abs',3600);
ylabel('Z Force (N)');
grid on;

subplot(212)
plot(tout,[rho,mf./Vf,mf]);
legend('Ambient Density','Profiler Density','Profiler Mass');
tticklabel('abs',3600);
ylabel('Density (kg/m^3); Mass (kg)');
grid on;

% Energy usage for thrusting.  Assuming Tecnadyne 560 operating mostly at 
% Bollard: From tecnadyne data: 
% thrust [N] = 1.4546*(power [W])^(2/3) + -7.4886
% @@@ this fit is for Bollard in the forward condition only.  
pwr = ((abs(Zthrust) + 7.4886*(Zthrust~=0))/1.4546).^(3/2);
warning('Using Tecnadyne thruster model for power estimates!');

figure(3); clf reset;
subplot(211)
plot(tout,pwr);
tticklabel('abs',3600);
ylabel('Thruster Power (W)');
grid on;

subplot(212)
J2KWH = 1/1000/3600;
plot(tout,[NaN; J2KWH*cumsum(pwr(2:end).*diff(tout))]);
tticklabel('abs',3600);
ylabel('Thruster Energy (kWh)');
grid on;



% For proposal.
figure(4); clf reset;
%set(gcf,'position',[985   303   935   554]);
set(gcf,'position',[1377         303         543         554]);
subplot(311)
plot(tout,yout(:,2));
tsupr = teout(ieout == strmatch('SUPRdriver',strvcat({prm.components.name})));
tsuprstart = tsupr(1:2:end);
tsuprstop = tsupr(2:2:end);
zsupr = NaN*tout;
for n=1:length(tsuprstart)
  ii = (tout >= tsuprstart(n) & tout < tsuprstop(n));
  zsupr(ii) = yout(ii,2);
end
line(tout,zsupr,'color','r');
set(gca,'ydir','reverse');
ylabel('Depth (m)');
legend('Ascent/Descent','Filtering','location','NW')
grid on;

subplot(312)
plot(tout,[rho,mf./Vf]);
legend('Ambient Density','Profiler Density');
ylabel('Density (kg/m^3)');
grid on;
ylim([1010 1110]);

subplot(313)
J2WH = 1/3600;
plot(tout,pwr,tout,[NaN; J2WH*cumsum(pwr(2:end).*diff(tout))]);
tticklabel('abs',3600*2);
ylabel(sprintf('Propulsion Power (W)\nEnergy (Wh)'));
legend('Power','Energy');
grid on;
xlabel('Time (hh:mm:ss)');

ch = getch(gcf,'axes');
set(ch(2:3),'xticklabel',[]);
vstretch(0);


% For proposal - smaller version.
figure(5); clf reset;
%set(gcf,'position',[985   303   935   554]);
set(gcf,'position',[1325         447         543         317]);
subplot(211)
plot(tout,yout(:,2));
line(tout,zsupr,'color','r');
set(gca,'ydir','reverse');
ylabel('Depth (m)');
legend('Ascent/Descent','Filtering')
grid on;

subplot(212)
plot(tout,pwr,tout,[NaN; J2WH*cumsum(pwr(2:end).*diff(tout))]);
tticklabel('abs',3600*2);
ylabel(sprintf('Propulsion Power (W)\nEnergy (Wh)'));
legend('Power','Energy');
grid on;
xlabel('Time (hh:mm:ss)');

ch = getch(gcf,'axes');
set(ch(2),'xticklabel',[]);
vstretch(0);

% For isopycnal RAFOS float study.
figure(6); clf reset;
subplot(211)
plot(tout,yout(:,2));
set(gca,'ydir','reverse');
ylabel('Depth (m)');
grid on;

subplot(212)
pden = sw_pden(S,theta-273.15,p*1e-4,0);
plot(tout,[rho,pden,mf./Vf]);
legend('Ambient Density','Ambient Pot. Density','Profiler Density');
ylabel('Density (kg/m^3)');
grid on;
ylim([1010 1110]);
tticklabel('abs',3600*1);
xlabel('Time (hh:mm:ss)');

