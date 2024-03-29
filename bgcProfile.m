function profile = bgcProfile(const)
% Create a profile compatible with simulation from data.
%
% Revision History
% 2012-12-27    mvj    Created placeholder function.
% 2013-01-04    mvj    Added Levitus 1982 World average profile.


% @@@ Placeholder profile.
% profile.z = [-1000 0 0:10:6000 11000];
% profile.theta = 273.15 + [25 25 linspace(25,2,length(profile.z)-3) 2]; % [K, in situ] 
% profile.rho = [0 0 linspace(1030,1050,length(profile.z)-3) 2000]; % [kg/m^3, in situ]
% profile.p = const.atm + [0 0.5*const.g ...
%       * cumsum((profile.rho(1:end-1)+profile.rho(2:end)).*(profile.z(2:end)-profile.z(1:end-1)))]; % [Pa]

% Levitus 1982 World Ocean
levitus82World.z = [ ...
0 ... 
10 ... 
20 ... 
30 ... 
50 ... 
75 ... 
100 ... 
125 ... 
150 ... 
200 ... 
250 ... 
300 ... 
400 ... 
500 ... 
600 ... 
700 ... 
800 ... 
900 ... 
1000 ... 
1100 ... 
1200 ... 
1300 ... 
1400 ... 
1500 ... 
1750 ... 
2000 ... 
2500 ... 
3000 ... 
3500 ... 
4000 ... 
4500 ... 
5000 ... 
5500 ...
]';

levitus82World.S = [ ...
34.63 ... 
34.69 ...
34.75 ... 
34.80 ... 
34.90 ... 
34.97 ... 
35.02 ... 
35.03 ... 
35.03 ... 
34.99 ... 
34.93 ... 
34.87 ... 
34.76 ... 
34.68 ... 
34.63 ... 
34.61 ... 
34.60 ... 
34.60 ... 
34.61 ... 
34.63 ... 
34.65 ... 
34.67 ... 
34.68 ... 
34.70 ... 
34.72 ... 
34.74 ... 
34.75 ... 
34.74 ... 
34.74 ... 
34.73 ... 
34.73 ... 
34.73 ... 
34.72 ...
]';

levitus82World.Theta = [ ...
18.10 ...
18.03 ...
17.89 ...
17.67 ...
17.00 ...
15.99 ...
14.99 ...
14.06 ...
13.24 ...
11.84 ...
10.75 ...
9.88 ...
8.46 ...
7.28 ...
6.38 ...
5.65 ...
5.04 ...
4.54 ...
4.12 ...
3.78 ...
3.50 ...
3.26 ...
3.04 ...
2.85 ...
2.47 ...
2.18 ...
1.79 ...
1.50 ...
1.27 ...
1.08 ...
0.91 ...
0.90 ...
1.02 ...
]';

levitus82World.db = sw_pres(levitus82World.z,0);
levitus82World.T = sw_temp(levitus82World.S,levitus82World.Theta,levitus82World.db,0);
levitus82World.rho = sw_dens(levitus82World.S,levitus82World.T,levitus82World.db);

% Select profile and add surface and hadal depth.
profile.z = [-1000; 0; levitus82World.z; 11000];
profile.theta = [25; 25; levitus82World.T; 2] + 273.15;
profile.rho = [0; 0; levitus82World.rho; sw_dens(levitus82World.S(end),levitus82World.T(end),11000)];
profile.p = const.atm + [0; 0; levitus82World.db*1e4; sw_pres(11000,0)*1e4];

