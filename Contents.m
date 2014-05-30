% Clio vehicle simulation.
%
% Revision History
% 2014-05-28    mvj    Added this documentation.
%
%
% Quick start:
% ------------
% At the matlab prompt, cd to the directory containing this file.
% Then either create a soft link (Linux) or create a copy of the 
% of the vehicle/simulation parameters file you would like to run.
% For example, to run the feasibility study vehicle:
% 
% >> !ln -s ./vehicle/20130200_proposal.m bgcParam.m  % Linux
% >> copyfile('./vehicle/20130200_proposal.m','./bgcParam.m')  % Windows
%
% then execute:
% >> bgc
%
% This will run a simulation of a Clio dive.  The vehicle 
% itself and all simulation parameters are defined in 
% bgcParam.m
%
% The only output as of 2015-05-28 is a series of plots.  
% No log is produced, but all variables are dumped to 
% the workspace and can be saved manually as a .mat file
% if desired.
%
%
% Function library:
% -----------------
%
% documentation pending...
%
% bgcAddComponent.m
% bgcBruntVaisala.m
% bgcBulkParam.m
% bgcEventBounds.m
% bgcEventControl.m
% bgcEventDropWeight.m
% bgcEventFilter.m
% bgcEventNone.m
% bgcEvents.m
% bgcFeedback.m
% bgcF.m
% bgcInitComponent.m
% bgcInSitu.m
% bgcIntegrator.m
% bgc.m
% bgcNeutral.m
% bgcParam.m
% bgcPlot.m
% bgcProfile.m
% bgcVolume.m
% seawater
% thrusters
% util
