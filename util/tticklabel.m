function tticklabel(varargin)
%TTICKLABEL  label x-axis/axes [HH:MM:SS]
%
%   TTICKLABEL('abs')
%   TTICKLABEL('abs',inc)
%   TTICKLABEL('rel',t0)
%   TTICKLABEL('rel',t0,inc)

% Revision History
% 2004          mvj    1.0    created
% 2005-07-12    mvj    1.0    'rel' does reasonable things for short intervals
% 2007-10-05    mvj    Modified xlabel to (hh:mm:ss)
% 2008-01-01    mvj    Reverts child order after running.


handle = get(gcf,'Children');
m = 0;
for ii = 1:length(handle)
  if strcmpi(get(handle(ii),'Type'),'axes') & ... 
	~strcmpi(get(handle(ii),'Tag'),'legend') & ...
	~strcmpi(get(handle(ii),'Tag'),'Colorbar')
    m = m+1;
    P = get(handle(ii),'Position');
    p(m) = P(2); % y-position
    ha(m) = handle(ii);
  end
end
M = m; % number of time subplots

[nul,ip] = sort(p); % sort ascending
ha = ha(ip([end:-1:1])); % sort descending 
if M == 2
  mm = 1:2;
else
  mm = M;
end

for m = 1:M
  axes(ha(m));

  switch varargin{1}
    case 'abs'
      if nargin > 1
	inc = varargin{2};
	XLIM = xlim;
	xtick = ceil(XLIM(1)/inc)*inc:inc:floor(XLIM(2)/inc)*inc;
	set(gca,'XTick',xtick);
      else
	xtick = get(gca,'XTick');
      end
      xtick = num2str(datestr(matlabtime(xtick), 13));
      set(gca,'XTickLabel',xtick);
      if any(m == mm)
	xlabel('Mission Day Time (hh:mm:ss)');	
      else
	xlabel('');    
      end
    case 'rel'
      T0 = varargin{2};
      if nargin > 2
	inc = varargin{3};
	XLIM = xlim;
	xtick = ceil(XLIM(1)/inc)*inc:inc:floor(XLIM(2)/inc)*inc;
	set(gca,'XTick',xtick);
      else
	xtick = get(gca,'XTick');
      end
      DT = diff(XLIM);
      if DT < 60 % [SS.FF]
	xtick = xtick-T0;
	set(gca,'XTickLabel',xtick);
	if any(m == mm)
	  xlabel('time [s]');
	else
	  xlabel('');
	end
      elseif DT < 3600 % [MM:SS]
	xtick = (xtick-T0)/60;
	set(gca,'XTickLabel',xtick);
	if any(m == mm)
	  xlabel('time [min]');
	else
	  xlabel('');
	end
      else % [HH:MM:SS]
	xtick = num2str(datestr(matlabtime(xtick-T0), 13));
	set(gca,'XTickLabel',xtick);
	if any(m == mm)
	  xlabel('Mission Time (hh:mm:ss)');
	else
	  xlabel('');
	end
      end
    otherwise
      error('Usage: tticklabel({''abs''|''rel'',T0},{inc})');
  end
  
end

set(gcf,'Children',handle);
