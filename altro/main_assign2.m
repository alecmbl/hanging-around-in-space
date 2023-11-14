%% Main script for the second assignment
clear
clc

%% Integration 

a = 10.9394e4;  % [km]
e = 0.3100;  % [-]
i = 46.2887;  OM = 0;  om = 0;  theta = 0;  % [rad]

k = 3;  m = 13;  % [-]

mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]

w_e = 15.04/3600;  % [deg/s]
thetag0 = 0; % [deg]
t0 = 0; % [s]

period = 2*pi*sqrt(a^3/mu_E); % Orbital period [1/s]
tspan = linspace( 0, 1*period, 10000 );

kep = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(theta), mu_E];
[r0, v0] = kep2car(kep);
y0 = [r0; v0];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[T, Y] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );


%% Ground Track

[~,~,lon,lat] = groundTrack(y0,thetag0,tspan,w_e,mu_E,t0);

wrap = wrapTo180(lon);

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap,lat,'g','linestyle','none','marker','.')
hold on
plot(wrap(1,1),lat(1,1),'*',wrap(end),lat(end),'*')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]'); ylabel('Latitude [deg]'); 
title('Ground Track');
legend('Ground Track', 'Start', 'End');
set(gca,'YDir','normal')

%% REPEATING GROUND TRACK

[a_rep,T_rep] = RepeatGroundTrack(k,m,w_e,mu_E);

t_rep = linspace( 0, 15*T_rep, 10000 );

[r0,v0] = kep2car([a_rep,e,i,OM,om,theta,mu_E]);
y0 = [r0; v0];

[alpha_rep,delta_rep,lon_rep,lat_rep] = groundTrack(y0,thetag0,t_rep,w_e,mu_E,t0);

wrap_rep = wrapTo180(lon_rep);

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap_rep,lat_rep,'r','linestyle','none','marker','.')
hold on
plot(wrap_rep(1,1),lat_rep(1,1),'o',wrap_rep(end),lat_rep(end),'o')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]'); ylabel('Latitude [deg]'); 
title('Repeating Ground Track');
legend('Repeating Ground Track', 'Start', 'End');
set(gca,'YDir','normal')




















