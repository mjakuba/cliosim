function [value,isterminal,direction] = bgcEventNone(t,y,varargin)
% For use with ODE45.  An event that will never be true.
%
% Revision History
% 2012-12-28    mvj    Created.
% 2019-04-09    mvj    Added varargin to allow passing arguments to other functions in bgcF.m


value = 1;
isterminal = 1;
direction = 0;

