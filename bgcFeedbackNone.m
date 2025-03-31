function [Zthrust,dizedt,diztedt] = bgcFeedbackPIV(t,zt,z,ize,izte,zg,prm)
% Revision History
% 2025-03-31    mvj    Created.

[Zthrust,dizedt,diztedt] = deal(0);
