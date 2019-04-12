function [Zthrust,dizedt,diztedt] = bgcFeedbackPIV(t,zt,z,ize,izte,zg,prm)
% Revision History
% 2019-04-09    mvj    Create.

% Parse out arguments.
Kp = prm{1};
Kv = prm{2};
Ki = prm{3};
ztmax = prm{4}; % this effectively sets desired transit speed.
Zmax = prm{5};

% Feed-forward terms.
Zballast = prm{6};
Zww = prm{7};

% velocity command
ze = zg-z;
ztg = Kp*ze;

% Saturate.   This will set maximum transit speed.
if abs(ztg) > ztmax
  ztg = sign(ztg)*ztmax;
end

% velocity error
zte = ztg - zt;

% Integral windup.  Output units are N.  This is operating on the velocity signal, 
% not the position signal.  Divisor is arbitrary.
% saturation implementation for external integration.  If above saturation and zte
% would further increase, set d/dt izte to 0.  Vice versa for below saturation.
dizedt = ze; % no saturation.  This term is not used.
Zsat = Zmax*0.5;
if Ki > 0
  if Ki*izte > Zsat && zte > 0
    diztedt = 0;
  elseif Ki*izte < -Zsat && zte < 0
    diztedt = 0;
  else
    diztedt = zte;  % moving away from saturation.
  end
end
%diztedt = zte; % @@@@ NEUTER saturation!

% Feedforward.
% @@@ Second term here is actually in the feedback path and seems to be a bad idea.
% But first term is needed at very start of sim to avoid issue with bounds (vehicle floats up before integral term
% does its work).
Zff = -Zballast;
%Zff = -Zballast + -Zww*ztg*abs(ztg);  % bad idea apparently!
%Zff = -Zballast + -Zww*zte*abs(zte);  %This might be ok, but unclear if it does anything useful
%[zt ztg -Zww*ztg*abs(ztg)]

% Feedback (PI on velocity).  I term is hopeless without FF.
Zthrust = Zff + Ki*izte + Kv*zte;


% Saturate.
if abs(Zthrust) > Zmax
  Zthrust = sign(Zthrust)*Zmax;
end


%[t/1000,zg/1000,z/1000,Zthrust,Zff,Ki*izte,Kv*zte,zt,ztg]
%[zg/1000,z,Zthrust,izte,ztg,zt]
