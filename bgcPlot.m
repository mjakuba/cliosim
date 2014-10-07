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
ie = ones(size(ieout));
% Controller is always last in component list, but actually we should add
% a component that is the filtering operation - an event function already 
% already exists bgcEventFilter - add a component that represents the
% filtering operation.
ie(ieout == length(prm.components)) = NaN; 
ie(1:2) = NaN;
zsupr = interp1(teout+0.1*rand(length(teout),1),ones(size(teout)).*ie,tout); % add jitter to accommodate interp1 uniqueness requirement.
zsupr = zsupr.*yout(:,2);
line(tout,zsupr,'color','r');
set(gca,'ydir','reverse');
ylabel('Depth (m)');
legend('Ascent/Descent','Filtering',4)
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
legend('Power','Energy',2);
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
legend('Ascent/Descent','Filtering',4)
grid on;

subplot(212)
plot(tout,pwr,tout,[NaN; J2WH*cumsum(pwr(2:end).*diff(tout))]);
tticklabel('abs',3600*2);
ylabel(sprintf('Propulsion Power (W)\nEnergy (Wh)'));
legend('Power','Energy',2);
grid on;
xlabel('Time (hh:mm:ss)');

ch = getch(gcf,'axes');
set(ch(2),'xticklabel',[]);
vstretch(0);
