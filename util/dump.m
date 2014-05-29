function dump
% DUMP  Dumps all variables in calling workspace to base workspace.
%
% Useful as an alternative to 'keyboard' while developing functions 
% that ultimately will return only a small portion of their scope.
%
% Revision History
% 2010-10-16    mvj    Created.


v = evalin('caller','who');
for n=1:length(v)
  assignin('base',v{n},evalin('caller',v{n}));
end
