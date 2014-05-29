function ch = getch(hf,varargin)
%ch = getch(hf)
%ch = getch(hf,'axes') Get only children of a specific type (e.g. axes)
%ch = getch(hf,'line','marker','.',...) Gets children of a specific type 
%    (e.g. lines) having some specific property (e.g. 'marker' set to '.')
%
% To get colorbars or legends, use
% ch = getch(hf,'axes','tag','Colorbar') for colorbars
% ch = getch(hf,'axes','tag','legend') for legends
%
% Revision History
% 2006-09-26    mvj    Created.
% 2007-12-11    mvj    Almost never use multi-types.  Changed additional
%                      args to specific properties, selection criteria that
%                      I use often.
% 2008-06-12    mvj    Added note about grabbing legends and colorbars.


ch0 = get(hf,'children');
ch = [];

if nargin > 1
  type = varargin{1};

  % When we grab for axes, usually we don't want legends and
  % colorbars.  Axes have empty tags.
  if strcmp(type,'axes') ...
	& isempty(strmatch('tag',varargin{2:2:end},'exact'))
    varargin{end+1} = 'tag';
    varargin{end+1} = '';
    nnargin = nargin+2;
  else
    nnargin = nargin;
  end
  
  for n = 1:length(ch0)
    
    % Check if child is of the desired type.
    if strmatch(get(ch0(n),'type'),type)
      ch = [ch; ch0(n)];
      
      % Check if child has the desired properties.  If not,
      % delete it from the list.
      if nnargin > 2
	for nn = 2:2:nnargin-1
	  prop0 = varargin{nn};
	  value0 = varargin{nn+1};
	  value = get(ch0(n),prop0);
	  if length(value) ~= length(value0) ...
		|| ~all(value == value0)
	    ch = ch(1:end-1);
	    break;
	  end
	end
	
      end
      
    end
  end

else
  ch = ch0;
end
