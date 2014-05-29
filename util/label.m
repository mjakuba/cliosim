function varargout = label(x,y,varargin)
%label(x,y)
%label(x,y,{z|title})
%label(x,y,z,title)
%h = label(...)
%
% Revision History
% 2007-10-26    mvj    Created.

h(1) = xlabel(x);
h(2) = ylabel(y);

if nargin == 3 
  if all(get(gca,'View') == [0 90])
    h(3) = title(varargin{1});
  else
    h(3) = zlabel(varargin{1});
  end
elseif nargin == 4
  h(3) = zlabel(varargin{1});
  h(4) = title(varargin{2});
end

if nargout == 1
  varargout{1} = h;
end
