function kep = car2kep(rr,vv,mu)
%
% Transformation from Cartesian state to orbital elements
% 
% PROTOTYPE:
%   kep = rv2kp(rr,vv,mu)
%
% DESCRIPTION:
%   Conversion function from Cartesian coordinates to Keplerian elements.
% 
% INPUT:
%   rr    [3x1] Position vector [km]
%   vv    [3x1] Velocity vector [km/s]
%   mu    [1x1] Gravitational parameter [km^3/s^2]
% 
% OUTPUT:
%   a     [1x1] Semi-major axis [km]
%   e     [1x1] Eccentricity [-]
%   i     [1x1] Inclination [rad]
%   Om    [1x1] RAAN [rad]
%   om    [1x1] Pericenter anomaly [rad]
%   theta [1x1] True anomaly [rad]
%
%   kep   [1x6] Vector that contains all Keplerian elements
%
%
%   Commented part at the end of the code gives the output of i,Om,om ...
%   and theta in [deg] instead of [rad]
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

% tollerance to check if the value is more/less than zero due to numerical
% errors that can modify the results
toll = 1e-5;

% cartesian versors
ii = [1; 0; 0];
jj = [0; 1; 0];
kk = [0; 0; 1];

% position and velocity vector norm
r = norm(rr);
v = norm(vv);

% semi-major axis [a]
a = 1/((2/r) - (v^2/mu));

% eccentricity [e]
ee = 1/mu * ( ( v^2 - (mu/r) ) * rr - dot( rr, vv ) * vv );
e = norm(ee);

% angular momentum [h]
hh = cross(rr,vv);
h = norm(hh);

%inclination [i]
i = acos(dot(hh,kk)/h);

% true anomaly [N]
nn = cross(kk,hh)/norm(cross(kk,hh));

% right ascension of the ascending node [Om]
if i < toll

     i = 0;
     nn = ii;
     Om = 0;

elseif dot(nn,jj) >= toll
       Om = real(acos(dot(ii,nn)));
    else
       Om = 2*pi - real(acos(dot(ii,nn)));
end
    
% argument of pericenter [om]
if e < toll

    e = 0;
    ee = nn;
    om = 0;

    % true anomaly [theta]
    if dot(vv,rr) >= -toll
        theta = acos(dot(rr,ee)/r);
    elseif dot(vv,rr) < -toll
        theta = 2*pi - acos(dot(rr,ee)/r);
    end

else if dot(ee,kk) > -toll
        om = acos(dot(nn,ee)/e);
     else 
        om = 2*pi - acos(dot(nn,ee)/e);
     end

     % true anomaly [theta]
     if dot(vv,rr) >= -toll
         theta = acos(dot(rr,ee)/(r*e));
     elseif dot(vv,rr) < -toll
         theta = 2*pi - acos(dot(rr,ee)/(r*e));
     end

end   

if theta < toll
    theta = 0;
end

kep = [a, e, i, Om, om, theta];

% results in [deg]
% i = i*180/pi;
% Om = Om*180/pi;
% om = om*180/pi;
% theta = theta*180/pi;

end



