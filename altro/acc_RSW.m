function a_RSW = acc_RSW(t,kep)
%
%
% DESCRIPTION
%  
%
% INPUT:
%
% OUTPUT:
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

% Constants

mu_E = astroConstants(13);
mu_moon = astroConstants(20);
J2 = astroConstants(9);
R_planet = astroConstants(23);


% J2 Perturbation

a=kep(1);   % semi-major axis 
e=kep(2);   % eccentricity 
i=kep(3);   % inclination 
OM=kep(4);  % right ascension of the ascending node 
om=kep(5);  % arguument of periapsis 
th=kep(6);  % true anomaly 

% Set cartesian initial data:
[r,~] = kep2car([a,e,i,OM,om,th,mu_E]); 

% Distance of the s/c from the planet
rnorm = norm(r);

a_J2 = -3/2*J2*mu_E*R_planet^2/(rnorm^4);
car2RSW_vector = [ 1-(3*(sin(i)^2)*(sin(th+om)^2));...
             (sin(i)^2)*sin(2*(th+om));...
             sin(2*i)*sin(th+om)];
a_J2_RSW = a_J2*car2RSW_vector;
    
% Moon Perturbation

% Time traformation from seconds to days
day = t/(60*60*24);

% Distance of the moon from the planet
[r_moon] = ephMoon(day);
r_moon_norm = norm(r_moon);

% Distance between s/c and moon
r_sc_moon = r_moon' - r;
r_sc_moon_norm = norm(r_sc_moon);

a_moon = mu_moon*(r_sc_moon/(r_sc_moon_norm^3) - r_moon'/(r_moon_norm^3));
% Rotation matrix from XYZ (cartesian) frame to RSW frame, 
car2RSW_matrix = [-sin(OM)*cos(i)*sin(om+th)+cos(OM)*cos(om+th), cos(OM)*cos(i)*sin(om+th)+sin(OM)*cos(om+th), sin(i)*sin(om+th);
             -sin(OM)*cos(i)*cos(om+th)-cos(OM)*sin(om+th), cos(OM)*cos(i)*cos(om+th)-sin(OM)*sin(om+th), sin(i)*cos(om+th);
                          sin(OM)*sin(i), -cos(OM)*sin(i), cos(i)];
a_Moon_RSW=car2RSW_matrix*a_moon;


% Total Perturbation
a_RSW = a_J2_RSW + a_Moon_RSW;















