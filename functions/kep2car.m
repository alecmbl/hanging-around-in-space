function [r, v] = kep2car(kep)
%
% Conversion from Keplerian elements to Cartesian state vectors
%
% PROTOTYPE:
%   [r, v] = kep2car(kep)
%
% DESCRIPTION:
%   Conversion function from Keplerian elements to Cartesian coordinates. Angles in radians.
%
% INPUT:
%   kep = [a, e, i, OM, om, th, mu] [1x7] Keplerian elements vector
%   a                               [1x1] Semi-major axis [km]
%   e                               [1x1] Eccentricity [-]
%   i                               [1x1] Inclination [rad]
%   OM                              [1x1] RAAN [rad]
%   om                              [1x1] Pericentre anomaly [rad]
%   th                              [1x1] True anomaly [rad]
%   mu                              [1x1] Gravitational parametre [km^3/s^2]
%
% OUTPUT:
%   r[3x1] Position vector [km]
%   v[3x1] Velocity vector [km/s]
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


%Definition of keplerian elements
a  = kep(1) ;            %[km]           Semi-major axis
e  = kep(2) ;            %[-]            Eccentricity
i  = kep(3) ;            %[rad]          Inclination
OM = kep(4) ;            %[rad]          RAAN
om = kep(5) ;            %[rad]          Pericentre anomaly
th = kep(6) ;            %[rad]          True anomaly
mu = kep(7) ;            %[km^3/s^2]     Gravitational parametre

% Calculation of semi-latus rectum and norm of position vector
p = a .* ( 1 - e.^2 );
R = p ./ ( 1 + e .* cos(th) );

% Definition of position and velocity vectors in perifocal reference
R_pf = R .* [ cos(th); sin(th); 0];
V_pf = sqrt( mu/p ) .* [ -sin(th); e + cos(th); 0];

% Perifocal to inertial reference rotation matrices
R_OM = [ cos(OM) sin(OM)    0  ; -sin(OM) cos(OM)   0   ;   0     0       1  ]; % RAAN rotation matrix
R_i  = [    1       0       0  ;     0    cos(i)  sin(i);   0  -sin(i) cos(i)]; % Pericentre anomaly rotation matrix
R_om = [ cos(om) sin(om)    0  ; -sin(om) cos(om)   0   ;   0     0       1  ]; % Inclination rotation matrix

T = R_OM' * R_i' * R_om';                                                       % Perifocal to inertial rotation matrix

% Calculation of position and velocity vectors in cartesian coordinates
r = T * R_pf;           % [km] Position vector in inertial reference frame
v = T * V_pf;           % [km/s] Velocity vector in inertial reference frame


