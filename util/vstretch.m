function FACTOR = vstretch(FACTOR,varargin)
% vstretch(FACTOR)
% vstretch(FACTOR,ha)
% FACTOR = vstretch(...) 
%
% Vertically stretch subplots to better fill
% a figure window.
%
% If figure contains more than one row of subplots,
% specify the plots to apply this to with the optional
% HA argument.
%
% Use FACTOR=0 to apply auto-stretch - good for subplots 
% that all have the same abscissa.  In this case the 
% FACTOR applied is returned.
%
% Note - behaves somewhat differently from hstretch in that
% does not preserve space between subplots.
%
% Revision History
% 2009-10-12    mvj    Created from hstretch.m
% 2010-08-03    mvj    Fixed for odd #'s.  Thought I had done that earlier...


if nargin < 2
  ch = getch(gcf,'axes');
else
  ch = varargin{1};
end

for n = 1:length(ch);
  p(n,:) = get(ch(n),'position');
end

% order subplots bottom to top
[nul,is] = sort(p(:,2));
ch = ch(is);
p = p(is,:);

if FACTOR==0.0
  dy = p(2,2) - (p(1,2)+p(1,4));  % empty y-space.
  FACTOR = dy/p(1,4);
end

N = length(ch);
if mod(N,2) % odd number of subplots.
  mm = ceil(N/2);
  set(ch(mm),'position',p(mm,:)+[0 -p(mm,4)*FACTOR/2 0 p(mm,4)*FACTOR]);
  for m = 1:mm-1
    set(ch(m),'position',p(m,:)+[0 -p(m,4)*FACTOR*(0.5) 0 p(m,4)*FACTOR]);
  end
  for m = mm+1:N
    set(ch(m),'position',p(m,:)+[0 p(m,4)*FACTOR*(0.5)-p(mm,4)*FACTOR 0 p(m,4)*FACTOR]);
  end
else % even number of subplots.
  mm = ceil(N/2);
  for m = 1:mm-1
    set(ch(m),'position',p(m,:)+[0 -p(m,4)*FACTOR*(0.5) 0 p(m,4)*FACTOR]);
  end
  for m = mm:N
    set(ch(m),'position',p(m,:)+[0 p(m,4)*FACTOR*(-0.5) 0 p(m,4)*FACTOR]);
  end
end


