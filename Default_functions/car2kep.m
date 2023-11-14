function kep = car2kep(rr,vv,mu)
%
% Transformation from Cartesian state to orbital elements
% 
% [a,e,i,Om, om, theta] = rv2kp(rr,vv,mu)
% 
% -------------------------------------------------------
% Input arguments:
% 
% rr         [3x1]    position vector             [km]
% vv         [3x1]    velocity vector             [km/s]
% mu         [1x1]    gravitational parameter     [km^3/s^2]
% 
% -------------------------------------------------------
% Output arguments:
% a          [1x1]    semi-major axis             [km]
% e          [1x1]    eccentricity                [-]
% i          [1x1]    inclination                 [rad]
% Om         [1x1]    RAAN                        [rad]
% om         [1x1]    pericenter anomaly          [rad]
% theta      [1x1]    true anomaly                [rad]
%
% Commented part at the end of the code gives the output of i,Om,om ...
% and theta in [deg] instead of [rad]
% 
% CONTRIBUTORS:
% Aditya Kumar
% Armelli Andrea
% Cambielli Alessandro
% Cappellari Giovanni
%
% Final version:January 2023
%

% cartesian versors
ii = [1 0 0];
jj = [0 1 0];
kk = [0 0 1];

% position and velocity vector norm
r = norm(rr);
v = norm(vv);

% semi-major axis  
a = 1/((2/r) - (v^2/mu));

% eccentricity
ee = 1/mu*((v^2 - (mu/r))*rr - dot(rr,vv)*vv);
e = norm(ee);

% angular momentum
hh = cross(rr,vv);
h = norm(hh);

% inclination
i = acos(dot(hh,kk)/h);

% node line
nn = cross(kk,hh)/norm(cross(kk,hh));

% if cicles used to check the range of the arcosine argument ...
% since arccosine isn't a bijective function

% right ascension of the ascending node 
if dot(nn,jj)>=0
    Om = acos(dot(ii,nn));
else
    Om = 2*pi - acos(dot(ii,nn));
end

% pericenter anomaly
if dot(ee,kk)>=0
    om = acos(dot(nn,ee)/e);
else
    om = 2*pi - acos(dot(nn,ee)/e);
end

% true anomaly
if dot(vv,rr)>=0
    theta = acos(dot(rr,ee)/(r*e));
else
    theta = 2*pi - acos(dot(rr,ee)/(r*e));
end

kep = [a,e,i,Om,om,theta];


% results in [deg]
% i = i*180/pi;
% Om = Om*180/pi;
% om = om*180/pi;
% theta = theta*180/pi;

end


