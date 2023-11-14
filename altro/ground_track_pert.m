function [alpha_pert,delta_pert,lon_pert,lat_pert] = ground_track_pert(a, e, i, OM, om, th0 , thG_t0, t_vect, omegaE, muP, t0)
%
% Groundtrack parameters computation
%
% PROTOTYPE:
%   [alpha,delta,lon,lat] = groundTrack(a, e, i, OM, om, th0 , thG_t0, t_vect, omegaE, muP, t0)
%
% DESCRIPTION:
%   Computes alpha, delta and longitude, latitude to plot the
%   ground track of a S/C orbiting around a planet of constant muP
%
% INPUT:
%   a                               [1x1]    Semi-major axis [km]
%   e                               [1x1]    Eccentricity [-]
%   i                               [1x1]    Inclination [rad]
%   OM                              [1x1]    RAAN [rad]
%   om                              [1x1]    Pericentre anomaly [rad]
%   th0                             [1x1]    True anomaly [rad]
%   thG_t0                          [1x1]    Earth's initial anomaly [rad]
%   t_vect                          [1xn]    Time vector [s]
%   omegaE                          [1x1]    Earth's angular velocity [rad/s]
%   muP                             [1x1]    Gravitational parametre [km^3/s^2]
%   t0                              [1x1]    Initial time
%
% OUTPUT:
%   alpha                           [1xn]   Right ascension [rad]
%   delta                           [1xn]   Declination [rad]
%   lon                             [1xn]   Longitude [rad]
%   lat                             [1xn]   Latitude [rad]
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
[ time, state ] = ode113( @(t,s) ode_2bp_perturbed(s,t), t_vect, [r0;v0], options);

% Parameters computation for ground track in time
for i = 1:length(time)
    r(i) = norm(state(i,1:3));                                      % [km]  position vector norm
  
    delta_pert(i) = asin(state(i,3)./r(i));                              % [rad] declination
    
    % If cycle to check for alpha sign
    if state(i,2)./r(i) > 0
        alpha_pert(i) = acos((state(i,1)./r(i))/cos(delta_pert(i)));          % [rad] right ascension
    else
        alpha_pert(i) = 2*pi - acos((state(i,1)./r(i))/cos(delta_pert(i)));   % [rad] right ascension
    end
    
    thG(i) = thG_t0 + omegaE.*(time(i)-t0);
    lon_pert(i) = alpha_pert(i) - thG(i);                                     % [rad] longitude
    
    lat_pert = delta_pert;                                                    % [rad] latitude
end