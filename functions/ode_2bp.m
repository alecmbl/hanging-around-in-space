function dy = ode_2bp( ~, y, mu )
%
% ode_2bp ODE system for the unperturbed two-body problem (Keplerian motion)
%
% PROTOTYPE:
%   dy = ode_2bp( ~, y, mu )
%
% DESCRIPTION:
%   Returns the ode function of the unperturbed two body problem, 
%   given as input the Time(can be omitted), the State Vector and 
%   the Gravitational Parameter of the primary
%
% INPUT:
%   t  [1x1] Time (can be omitted, as the system is autonomous) [s]
%   y  [6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [km, km/s]
%   mu [1x1] Gravitational parameter of the primary [km^3/s^2]
%
% OUTPUT:
%   dy [6x1] Derivative of the state [km/s^2, km/s^3]
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

% Position and velocity of the body
r = y(1:3);
v = y(4:6);

% Distance from the primary body
rnorm = norm(r);

% Set the derivatives of the state
dy = [ v ; (-mu/rnorm^3)*r ];

end