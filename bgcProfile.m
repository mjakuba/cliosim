function profile = bgcProfile(const,varargin)
% Create a profile compatible with simulation from data.
%
% prm.const = bgcConst
% profile = bgcProfile(prm.const)  % default Levitus 1982 world average
%         = bgcProfile(prm.const,'isopycnal/isothermal')  % S varies with depth
%         = bgcProfile(prm.const,'isopycnal/isohaline')   % T varies with depth
%         = ...
%         = bgcProfile(prm.const,'user',S_PSU,T_degC,P_db)
%
% Revision History
% 2012-12-27    mvj    Created placeholder function.
% 2013-01-04    mvj    Added Levitus 1982 World average profile.
% 2022-03-12    mvj    Added some archetypal profiles and the option for an externally specified one.


if nargin == 2
    profile_type = varargin{1};
elseif nargin == 4
    profile_type = 'user';
    user.S = varargin{1}(:);
    user.T = varargin{2}(:);
    user.db = varargin{3}(:);
    in = isnan(user.S) | isnan(user.T) | isnan(user.db);
    user.S = user.S(~in);
    user.T = user.T(~in);
    user.db = user.db(~in);
    assert(all(user.db >= 0),'User specified profile cannot have negative pressures.');
else
    profile_type = 'levitus82World';
end

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

levitus82World.db = sw_pres(levitus82World.z,const.lat);
levitus82World.T = sw_temp(levitus82World.S,levitus82World.Theta,levitus82World.db,0);
levitus82World.rho = sw_dens(levitus82World.S,levitus82World.T,levitus82World.db);

% compose isopycnal profile.  Also isohaline.
if strcmp(profile_type,'isopycnal/isohaline')
    isopycnal.z = [0:10:11000]';
    pden0 = 1030;
    p0 = 0;
    S0 = 35;
    fprintf(1,'Generating isopycnal/isohaline profile...');
    for n = 1:length(isopycnal.z)
        z = isopycnal.z(n);
        p = sw_pres(z,const.lat);
        f = @(temp,S,p,p0,pden0) abs(sw_pden(S,temp,p,p0)-pden0);
        isopycnal.T(n,1) = fminsearch(@(temp) f(temp,S0,p,p0,pden0),10);
        isopycnal.db(n,1) = p;
    end
    isopycnal.S = S0*ones(size(isopycnal.z));
    isopycnal.rho = sw_dens(isopycnal.S,isopycnal.T,isopycnal.db);
    isopycnal.sigma0 = sw_pden(isopycnal.S,isopycnal.T,isopycnal.db,0);
    assert(range(isopycnal.sigma0) < 1e-4);
    fprintf(1,'Done.');
end

% compose isopycnal profile.  Also isothermal (in situ)
if strcmp(profile_type,'isopycnal/isothermal')
    isopycnal.z = [0:10:11000]';
    pden0 = 1030;
    p0 = 0;
    T0 = 10;
    fprintf(1,'Generating isopycnal/isothermal profile...');
    for n = 1:length(isopycnal.z)
        z = isopycnal.z(n);
        p = sw_pres(z,const.lat);
        f = @(temp,S,p,p0,pden0) abs(sw_pden(S,temp,p,p0)-pden0);
        isopycnal.S(n,1) = fminsearch(@(salt) f(T0,salt,p,p0,pden0),35);
        isopycnal.db(n,1) = p;
    end
    isopycnal.T = T0*ones(size(isopycnal.z));
    isopycnal.rho = sw_dens(isopycnal.S,isopycnal.T,isopycnal.db);
    isopycnal.sigma0 = sw_pden(isopycnal.S,isopycnal.T,isopycnal.db,0);
    assert(range(isopycnal.sigma0) < 1e-4);
    fprintf(1,'Done.');
end


% Select profile and add surface and hadal depth.
switch profile_type
  case 'levitus82World'
    profile.z = [-1000; 0; levitus82World.z; 11000];
    profile.theta = [25; 25; levitus82World.T; 2] + 273.15;
    profile.S = [0; 0; levitus82World.S; levitus82World.S(end)];
    profile.rho = [0; 0; levitus82World.rho; sw_dens(levitus82World.S(end),levitus82World.T(end),11000)];
    profile.p = const.atm + [0; 0; levitus82World.db*1e4; sw_pres(11000,const.lat)*1e4];
  case {'isopycnal/isohaline','isopycnal/isothermal'}
    profile.z = [-100; 0; isopycnal.z];
    profile.theta = [25; 25; isopycnal.T] + 273.15;
    profile.S = [0; 0; isopycnal.S];
    profile.rho = [0; 0; isopycnal.rho];
    profile.p = const.atm + [0; 0; isopycnal.db*1e4];
    profile.sigma0 = [0; 0; isopycnal.sigma0];
  case 'pycnoclinic/isohaline'  % easiest to understand because there is no salinity effect.
    z = linspace(0,1000,300)';
    db = sw_pres(z,const.lat);
    S = 35*ones(size(z));
    T = linspace(25,2,300)';
    rho = sw_dens(S,T,db);
    assert(all(diff(sw_pden(S,T,db,0)) > 0));  % profile is stable.
    profile.z = [-100; 0; z];
    profile.S = [0; 0; S];
    profile.theta = [25; 25; T] + 273.15;
    profile.p = const.atm + [0; 0; db*1e4];
    profile.rho = [0; 0; rho];
  case 'pycnoclinic/positive-upward haloclinic'  % worst case for an isopycnal float
    z = linspace(0,1000,300)';
    db = sw_pres(z,const.lat);
    S = linspace(37,35,300)';
    T = linspace(25,2,300)';
    rho = sw_dens(S,T,db);
    assert(all(diff(sw_pden(S,T,db,0)) > 0));  % profile is stable.
    profile.z = [-100; 0; z];
    profile.S = [0; 0; S];
    profile.theta = [25; 25; T] + 273.15;
    profile.p = const.atm + [0; 0; db*1e4];
    profile.rho = [0; 0; rho];
  case 'pycnoclinic/positive-downward haloclinic'
    z = linspace(0,1000,300)';
    db = sw_pres(z,const.lat);
    S = linspace(35,37,300)';
    T = linspace(25,2,300)';
    rho = sw_dens(S,T,db);
    assert(all(diff(sw_pden(S,T,db,0)) > 0));  % profile is stable.
    profile.z = [-100; 0; z];
    profile.S = [0; 0; S];
    profile.theta = [25; 25; T] + 273.15;
    profile.p = const.atm + [0; 0; db*1e4];
    profile.rho = [0; 0; rho];
  case 'user'
    profile.z = [ -100; 0; sw_dpth(user.db,const.lat)];
    profile.theta = [25; 25; user.T] + 273.15;
    profile.S = [0; 0; user.S];
    profile.rho = [0; 0; sw_dens(user.S,user.T,user.db)];
    profile.p = const.atm + [0; 0; user.db*1e4];
  otherwise
    error('Unknown profile.');
end
 

