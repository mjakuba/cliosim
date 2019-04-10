function [Zthrust] = bgcFeedbackPID(t,zt,z,zg,prm)
% Revision History
% 2013-01-02    mvj    Created.
% 2019-04-08    mvj    rename to PID to allow for other control schemes.


% Parse out arguments.
Kp = prm{1};
Kd = prm{2};
Ki = prm{3};
Zff = prm{4};
Zmax = prm{5};

% Rate error.
zte = 0-zt;

% Position error.
ze = zg-z;

% Integral error.
ize = bgcIntegrator(t,ze);

% Feedback.
Zthrust = Kp*ze + Kd*zte + Ki*ize + Zff;

if abs(Zthrust) > Zmax
  Zthrust = sign(Zthrust)*Zmax;
end
