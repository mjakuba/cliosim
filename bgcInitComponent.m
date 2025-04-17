function component = bgcInitComponent(name)
% Revision History
% 2012-12-31    mvj    Created.


component.name = name;
component.m = 0;
component.V = 0;
component.rho = 0;
component.alpha = 0;
component.chi = 0;
component.cp = 0;
component.active = false;
component.activate_time = 0; % simulation start for components that start active.
component.discharge_rate = 0; % [m^3/s]
component.eos = @bgcVolumeNone;
component.eventf = @bgcEventNone;
component.event_prm = {};
