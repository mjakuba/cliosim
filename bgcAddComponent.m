function components = bgcAddComponent(component,varargin)
% Add a component to the components structure.  This just makes sure the
% fields are in the right order and that none are missing.
%
% Revision History
% 2012-12-28    mvj    Created.

c0 = bgcInitComponent('');

fn = fieldnames(c0);
for f = 1:length(fn)
  if isfield(component,fn{f})
    c0.(fn{f}) = component.(fn{f});
  end
end

if nargin < 2
  components = c0;
else
  components = varargin{1};
  components(end+1) = c0;
end
  
