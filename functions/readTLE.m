function [kep, EpochY, EpochD] = readTLE(file, catalog)
%
% READTLE Read satellite ephemeris data from a NORAD two-line element (TLE) file.
% 
% PROTOTYPE:
%   [kep, EpochY, EpochD] = readTLE(file, catalog)
%
% INPUTS:
%   file    - Path to any standard two-line element file.
%   catalog - Optional array of NORAD catalog numbers for the satellites of
%             interest. The default action is to display data from every
%             satellite in the file.
% OUTPUTS:
%   kep: keplerian elements matrix
%       a  [1x1] Semi-major axis [km]
%       e  [1x1] Eccentricity [-]
%       i  [1x1] Inclination [rad]
%       OM [1x1] RAAN [rad]
%       om [1x1] Pericentre anomaly [rad]
%       th [1x1] True anomaly [rad]
%   EpochY [1x1] last two digits of the year [-]
%   Epochd [1x1] day of the year and fraction of the day [-]
%
% CONTRIBUTORS:
%   Aditya Kumar
%   Armelli Andrea
%   Cambielli Alessandro
%   Cappellari Giovanni
%
% Final version:January 2023
%
% -------------------------------------------------------------------------


if nargin < 2
    catalog = [];
  end

  fd = fopen(file,'r');
  if fd < 0, fd = fopen([file '.tle'],'r'); end
  assert(fd > 0,['Can''t open file ' file ' for reading.'])
    
  kep = [];
  EpochY = [];
  EpochD = [];

  q = 0;
  A1 = fgetl(fd);
  A2 = fgetl(fd);

  while ischar(A2)
    q = q + 1;
    satnum = str2num(A1(3:7));
    if isempty(catalog) || ismember(satnum, catalog)
      %Epoch
      EpochYY = str2num(A1(19:20));
      EpochY = [EpochY; EpochYY];
      EpochDD = str2num(A1(21:32));
      EpochD = [EpochD; EpochDD];

      %Inclination
      i = deg2rad(str2num(A2(9:16)));

      % RAANS
      OM = deg2rad(str2num(A2(18:25)));

      %Eccentricity
      e = str2num(['.' A2(27:33)]);

      %Argument of perigee
      om = deg2rad(str2num(A2(35:42)));

      %Mean Anomaly
      M = deg2rad(str2num(A2(44:51)));

      %Mean motion
      n = str2num(A2(53:63));

      %Semi-major axis
      T = 86400/n;
      a = ((T/(2*pi))^2*398.6e12)^(1/3)*1e-3;

      %True Anomaly
      options = optimset('Display','off', 'TolFun', 1e-13);
      fun = @(E) E -i*sin(E) -M;
      E = fzero(fun,pi,options);
      th = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
      th = wrapTo2Pi(th);
      
      kepp = [a, e, i, OM, om, th];
      kep = [kep; kepp];
      
    end

    A1 = fgetl(fd);
    A2 = fgetl(fd);
  end

  fclose(fd);
end
