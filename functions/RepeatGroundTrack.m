function [a_rep,T_rep] = RepeatGroundTrack(k,m,w_e,mu)
%
% Compute the semi-major axis for a repeating orbit.
%
% PROTOTYPE:
%   [a_rep,T_rep] = RepeatGroundTrack(k,m,w_e,mu)
%
% INPUT:
%   k     [1x1] Number of satellite revolutions [-]
%   m     [1x1] Number of earth revolutions [-]
%   mu    [1x1] Gravitational parameter of the primary [m^3/s^2]
%   w_e   [1x1] Earth's rotation velocity (eastwards) [deg/s]
%
% OUTPUT:
%   a_rep [1x1] Semi-major axis of the repeating orbit [km]
%   T_rep [1x1] Orbital period of the repeating orbit [s]
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

% set the new orbital period considering the number of satellite and
% earth's revolutions after which the groundtrack must repeat itself:
T_rep = 360/w_e*m/k;  %[s]

% derive the mean orbital motion n:
n = 2*pi/T_rep;  %[1/s]

% knowing the relation between the mean motion and the semi-major axis,
% it is possible to derive a:
a_rep = (mu/(n^2))^(1/3);  %[km]
