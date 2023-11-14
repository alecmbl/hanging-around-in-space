function dy = ode_2bp_perturbed(y, time)
%
% ode_2bp_perturbed ODE system for the perturbed two-body problem (Keplerian motion)
%
% PROTOTYPE:
%   dy = ode_2bp_perturbed(t, y)
%
% DESCRIPTION
%   Returns the ode function of the perturbed (J2 and moon perturbations)
%   two body problem, given as input the Time and the State Vector.
%
% INPUT:
%   t[1] Time  [s]
%   y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [km, km/s ]
%
% OUTPUT:
%   dy[6x1] Derivative of the state [km/s^2, km/s^3]
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

%% Constants

mu_E = astroConstants(13);
mu_moon = astroConstants(20);
J2 = astroConstants(9);
R_planet = astroConstants(23);

%% J2 Perturbation

% Position and velocity of the s/c wrt planet
r = y(1:3);
v = y(4:6);

% Distance of the s/c from the planet
rnorm = norm(r);

a_J2 = 3/2 * J2 * mu_E * R_planet^2 / rnorm^4 .*...
    [r(1) / rnorm*( 5 *(r(3) / rnorm)^2 - 1);...
    r(2) / rnorm*( 5 *(r(3) / rnorm)^2 - 1);...
    r(3) / rnorm*( 5 *(r(3) / rnorm)^2 - 3)];

%% Moon Perturbation

% Time traformation from seconds to days
day = time/(60*60*24);

% Distance of the moon from the planet
[r_moon] = ephMoon(day)';
% r_moon = r_moon';
r_moon_norm = norm(r_moon);

% Distance between s/c and moon
r_sc_moon = r_moon - r;
r_sc_moon_norm = norm(r_sc_moon);

a_moon = mu_moon * ( r_sc_moon ./ (r_sc_moon_norm^3) - r_moon ./ (r_moon_norm^3));

%% Total Perturbation
a = a_J2 + a_moon;

% Set the derivatives of the state
dy = [ v ; (-mu_E/rnorm^3).*r + a ];















