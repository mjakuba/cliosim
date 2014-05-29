function [value,isterminal,direction] = bgcEvents(t,y,prm)
% For use with ODE45.  Aggregates event functions associated with
% individual components.
%
% Revision History
% 2012-12-28    mvj    Created.

value = [];
isterminal = [];
direction = [];
for c = 1:length(prm.components)
  [v,ist,d] = prm.components(c).eventf(t,y,prm.components(c).event_prm{:});
  value = [value; v];
  isterminal = [isterminal; ist];
  direction = [direction; d];
end


