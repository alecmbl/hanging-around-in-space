function [alpha,delta,lon,lat] = ground_track_perturbed(a, e, i, OM, om, th0 , thG_t0, t_vect, omegaE, muP, t0)
%
% Groundtrack parameters computation for the perturbed problem.
%
% PROTOTYPE:
%   [alpha,delta,lon,lat] = groundTrack(a, e, i, OM, om, th0 , thG_t0, t_vect, omegaE, muP, t0)
%
% DESCRIPTION:
%   Computes alpha, delta and longitude, latitude to plot the
%   ground track of a S/C orbiting around a planet of constant muP
%
% INPUT:
%   a      [1x1]    Semi-major axis [km]
%   e      [1x1]    Eccentricity [-]
%   i      [1x1]    Inclination [rad]
%   OM     [1x1]    RAAN [rad]
%   om     [1x1]    Pericentre anomaly [rad]
%   th0    [1x1]    True anomaly [rad]
%   thG_t0 [1x1]    Earth's initial anomaly [rad]
%   t_vect [1xn]    Time vector [s]
%   omegaE [1x1]    Earth's angular velocity [rad/s]
%   muP    [1x1]    Gravitational parameter [km^3/s^2]
%   t0     [1x1]    Initial time
%
% OUTPUT:
%   alpha  [1xn]   Right ascension [rad]
%   delta  [1xn]   Declination [rad]
%   lon    [1xn]   Longitude [rad]
%   lat    [1xn]   Latitude [rad]
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


% Conversion from keplerian elements to cartesian representation for the
% initial state
kep = [a e i OM om th0 muP];
[r0, v0] = kep2car(kep);

% Set options for ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-13);

% Integration of the orbit
[ time, state ] = ode113( @(t,s) ode_2bp_perturbed(t,s), t_vect, [r0;v0], options);

% Parameters computation for ground track in time
for i = 1:length(time)
    r(i) = norm(state(i,1:3));                                      % [km]  position vector norm
  
    delta(i) = asin(state(i,3)./r(i));                              % [rad] declination
    
    % If cycle to check for alpha sign
    if state(i,2)./r(i) > 0
        alpha(i) = acos((state(i,1)./r(i))/cos(delta(i)));          % [rad] right ascension
    else
        alpha(i) = 2*pi - acos((state(i,1)./r(i))/cos(delta(i)));   % [rad] right ascension
    end
    
    thG(i) = thG_t0 + omegaE.*(time(i)-t0);
    lon(i) = alpha(i) - thG(i);                                     % [rad] longitude
    
    lat = delta;                                                    % [rad] latitude
end