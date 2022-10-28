load ~/veh/rafos_compressee/sensitivity/simOutput.mat

figure(1); clf reset;
for nn=1:length(simOutput)

    s = simOutput{nn};
    
    plot(s.tout,s.yout(:,2));
    hold on;

end
set(gca,'ydir','reverse');
plot(xlim,[400 400],'--');
plot(xlim,[350 350],'--');
ylabel('Depth (m)');
tticklabel('abs',3600*1);
xlabel('Time (hh:mm:ss)');
grid on;


figure(2); clf reset;
% using initial condition for target pot density; not quite right but ballast depth not archived in sim.
pden0 = sw_pden(s.S(1),s.theta(1)-273.15,s.p(1)*1e-4,0);
for nn=1:length(simOutput)

    s = simOutput{nn};

    pden = sw_pden(s.S,s.theta-273.15,s.p*1e-4,0);  % pden along float trajectory
    plot(s.tout,pden-pden0);
    hold on;

end
% this is the ambient potential density versus the target at the target depth.
ylabel('Potential Density Variation (kg/m^3)');
grid on;
tticklabel('abs',3600*1);
xlabel('Time (hh:mm:ss)');


% will probably want a version of this where the density difference is plotted too, that is, the difference in potential density and the float?

% Swift et al., 1994 is probably the way we want to analyse this.
% The method therein uses a Taylor Series expansion to estimate the difference between the target iso-whatever and the equilibrium iso-whatever as a function of water column stratification (this can be converted into depth given a profile).  This allows assessment of how the deviation varies spatially throughout an ocean.  They also assess how ballasting and variability in actual compressibility affect target depth vs. realized depth.  What I want to know is how variability in compressibility affects deviation in realized iso-whatever vs. target iso-whatever, as the float moves around.  Every float is ballasted individually, and this sim as written assumes that was done perfectly to hit the target isobar at the deployment site (which is determined from a target potential density and a profile).

% get target potential densities
