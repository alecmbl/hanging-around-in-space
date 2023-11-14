function [a_repJ2,T_repJ2] = RepeatGroundTrack_J2(k,m,w_e,mu_E,e,i,a_rep)
% 
% Function that computes the semi-major axis required to obtain a repeating
% ground track with the given values of k,m, taking into account the
% secular effect of J2
% 
% PROTOTYPE:
%   [a_repJ2,T_repJ2] = RepeatGroundTrack_J2(k,m,w_e,mu_E,e,i,a_rep) 
%
% INPUT:
%  k       [1x1] Number of revolutions of the s/c to obtain the ground track repetition [-]
%  m       [1x1] Number of revolutions of the planet to obtain the ground track repetition [-]
%  mu      [1x1] Gravitational parameter of primary body [km^3/s^2]
%  w_e     [1x1] Earth's spin rate [rad/s]
%  e       [1x1] S/c orbit's eccentricity [-]
%  i       [1x1] S/c orbit's inclination [rad]
%  a_rep   [1x1] Semi-major axis for the repetition of the unperturbed ground track [km]
%
% OUTPUT: 
%  a_repJ2 [1x1] Semi-major axis of the repeating orbit [km]
%  T_repJ2 [1x1] Orbital period of the repeating orbit [s]
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

R_e = astroConstants(23);  % mean equatorial Earth radius [km]
J2 = astroConstants(9);  % J2 parameter for perturbation

c = (-3/2)*((sqrt(mu_E)*J2*(R_e^2))/((1-e^2)^2)); % constant coefficient

% a_rep is used as initial guess for fzero
a_repJ2 = fzero(@(a) (m/k) - ((w_e - (c/a^(7/2))*cos(i))/...
    ((sqrt(mu_E/a^3))+(c/a^(7/2))*((5/2)*(sin(i))^2-2)-...
    ((c*sqrt(1-e^2))/a^(7/2))*(1-(3/2)*(sin(i))^2))),a_rep);

T_repJ2 = 2*pi*sqrt((a_repJ2^3)/mu_E);


















