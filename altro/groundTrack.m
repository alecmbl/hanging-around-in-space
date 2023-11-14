function [alpha,delta,lon,lat] = groundTrack(y0,thetag0,t,w_e,mu_E,t0)
% Compute the ground track of the orbit

% Set options for the ODE solver

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration

[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), t, y0, options );

thetag = thetag0 + w_e*(t - t0); % [deg]

x = zeros(length(Y),1);
y = zeros(length(Y),1);
z = zeros(length(Y),1);
r = zeros(length(Y),1);
delta = zeros(1,length(Y)); 
alpha = zeros(1,length(Y)); 

for i = 1:length(Y)
    x(i,1) = Y(i,1);
    y(i,1) = Y(i,2);
    z(i,1) = Y(i,3);
    r(i,1) = sqrt((x(i,1)^2) + (y(i,1)^2) + (z(i,1)^2));
    delta(1,i) = asind(z(i,1)/r(i,1)); % [deg]
    alpha(1,i) = atan2(y(i,1),x(i,1))*180/pi; % [deg]
end

lon = alpha - thetag; % [deg]
lat = delta; % [deg]




























end