# Clio Vehicle Simulation

A Matlab library for simulating Clio dives with a particular focus on compressibility effects.  See Contents.m for further details.

## Quick Start

At the Matlab prompt, cd to the directory containing this file.
Then either create a soft link (Linux) or create a copy of the 
of the vehicle/simulation parameters file you would like to run.
For example, to run the feasibility study vehicle:

Linux:
```
>> !ln -s ./vehicle/20130200_proposal.m bgcParam.m  # Linux
```

Windows:
```
>> copyfile('./vehicle/20130200_proposal.m','./bgcParam.m')  # Windows
```

Now run the simulation

```
>> bgc
```

This will run a simulation of a Clio dive.  The vehicle 
itself and all simulation parameters are defined in 
bgcParam.m

The only output as of 2015-05-28 is a series of plots.  
No log is produced, but all variables are dumped to 
the workspace and can be saved manually as a .mat file
if desired.
