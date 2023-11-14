function [dkep] = GaussPert (kep,mu_earth,t)
%
% GaussPert.m
% 
% PROTOTYPE:
%   [dkep] = GaussPert (kep,mu_earth,t)
%
% DESCRIPTION:
% 	This function implements the Gauss equation for the computation of 
% 	perturbed orbital parameters 
%
% INPUT:
%   kep      [6x1] Keplerian parameters vector
%	mu_earth [1x1] Earth gravitational constant [km^3/s^2]
%   t        [nx1] Time interval of integration [s]
%
% OUTPUT:
%   [dkep]   [6x1] Perturbed keplerian parameters
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

a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
th = kep(6);

b = a*sqrt(1-e^2);
p = b^2/a;
n = sqrt(mu_earth/a^3);
h = n*a*b;
r = p/(1+e*cos(th));
v = sqrt((2*mu_earth)/r-mu_earth/a);


a_RSW = a_car2RSW(t,kep);
a_R= a_RSW(1);
a_S= a_RSW(2);
a_W= a_RSW(3);

da = (2*(a^2)/h) * (e*sin(th)*a_R + (p/r)*a_S);
de = (1/h) * ( p*sin(th)*a_R + ( (p+r)*cos(th) + r*e )*a_S );
di = ( r*cos(th+om)/h ) * a_W;
dOM = ( r*sin(th+om)/(h*sin(i)) ) * a_W;
dom = (1/(h*e)) * ( -p*cos(th)*a_R + (p+r)*sin(th)*a_S ) - ( (r*sin(th+om)*cos(i)) / (h*sin(i)) )*a_W ;
dth =  (h/(r^2)) + (1/(e*h)) * ( p*cos(th)*a_R - (p+r)*sin(th)*a_S );

dkep = [da,de,di,dOM,dom,dth]';








