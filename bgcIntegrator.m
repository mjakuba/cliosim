function [izeout,tizeout] = bgcIntegrator(t,ze,varargin);

error('This function is deprecated in favor of letting ode45 handle all integration.');

persistent ize tize;

% Hard set integrator if requested.  [] resets.
if nargin > 2
  ize = varargin{1};
  tize = varargin{2};
end

% Integral error.
if ~isempty(tize)
  ize = sum([ize ze*(t-tize)]);
else
  ize = 0;
end
tize = t;

izeout = ize;
tizeout = tize;
