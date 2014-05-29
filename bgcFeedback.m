function [Zthrust] = bgcFeedback(t,zt,z,zg,Kp,Kd,Ki,Zff)
% Revision History
% 2013-01-02    mvj    Created.

% Rate error.
zte = 0-zt;

% Position error.
ze = zg-z;

% Integral error.
ize = bgcIntegrator(t,ze);

% Feedback.
Zthrust = Kp*ze + Kd*zte + Ki*ize + Zff;

