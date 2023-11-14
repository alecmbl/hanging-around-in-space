%%
clear,clc,close all;
addpath('time\');

% Nominal orbit:
a = 10.9394e4;  % semi-major axis [km]
e = 0.3100;  % eccentricity
i = 46.2887;  % inclination [deg]
OM = 0;  % right ascension of ascending node [deg]
om = 0;  % argument of perigee [deg]
th = 0;  % initial true anomaly [deg]

% Constants:
mu_E = astroConstants(13);  % gravitational parameter of the Earth [km^3/s^2]
mu_M = astroConstants(20);  % gravitational parameter of the Moon [km^3/s^2]
R_E = astroConstants(23);  % mean equatorial Earth radius [km]
J2 = astroConstants(9);  % J2 parameter for perturbation

T = 2*pi*sqrt((a^3)/mu_E);  % orbital period [s]

%% Nominal Ground Tracks

N = 100000;  % number of points in the timespan 
oneday = 24*60*60;
tendays = 10*oneday;

tspan1 = linspace(0, T, N);
tspan2 = linspace(0, oneday, N);
tspan3 = linspace(0,tendays, N);

t0 = 0;  % initial time [s]
thG_t0 = 0;  % Earth's initial anomaly [rad] 
omegaE = deg2rad(15.04)/3600;  % Earth's angular velocity [rad/s]

% Unperturbed Groundtrack
[alpha1,delta1,lon1,lat1] = ground_track(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan1, omegaE, mu_E, t0);
[alpha2,delta2,lon2,lat2] = ground_track(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan2, omegaE, mu_E, t0);
[alpha3,delta3,lon3,lat3] = ground_track(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan3, omegaE, mu_E, t0);

lon1 = rad2deg(lon1);
lat1 = rad2deg(lat1);
lon2 = rad2deg(lon2);
lat2 = rad2deg(lat2);
lon3 = rad2deg(lon3);
lat3 = rad2deg(lat3);
wrap1 = wrapTo180(lon1);
wrap2 = wrapTo180(lon2);
wrap3 = wrapTo180(lon3);

% Plot of unperturbed ground tracks
figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap1,lat1,'g','linestyle','none','marker','.');
plot(wrap1(1,1),lat1(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap1(end),lat1(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - one period');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Unperturbed','Start','End');
set(gca,'YDir','normal')

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap2,lat2,'g','linestyle','none','marker','.');
plot(wrap2(1,1),lat2(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap2(end),lat2(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - one day');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Unperturbed','Start','End');
set(gca,'YDir','normal')

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap3,lat3,'g','linestyle','none','marker','.');
plot(wrap3(1,1),lat3(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap3(end),lat3(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - ten days');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Unperturbed','Start','End');
set(gca,'YDir','normal')

% Perturbed Groundtrack
[alpha1_pert,delta1_pert,lon1_pert,lat1_pert] = ground_track_pert(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan1, omegaE, mu_E, t0);
[alpha2_pert,delta2_pert,lon2_pert,lat2_pert] = ground_track_pert(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan2, omegaE, mu_E, t0);
[alpha3_pert,delta3_pert,lon3_pert,lat3_pert] = ground_track_pert(a, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan3, omegaE, mu_E, t0);

lon1_pert = rad2deg(lon1_pert);
lat1_pert = rad2deg(lat1_pert);
lon2_pert = rad2deg(lon2_pert);
lat2_pert = rad2deg(lat2_pert);
lon3_pert = rad2deg(lon3_pert);
lat3_pert = rad2deg(lat3_pert);
wrap1_pert = wrapTo180(lon1_pert);
wrap2_pert = wrapTo180(lon2_pert);
wrap3_pert = wrapTo180(lon3_pert);

% Plot of unperturbed ground tracks
figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap1_pert,lat1_pert,'r','linestyle','none','marker','.');
plot(wrap1_pert(1,1),lat1_pert(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap1_pert(end),lat1_pert(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - one period');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Unperturbed','Start','End');
set(gca,'YDir','normal')

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap2_pert,lat2_pert,'r','linestyle','none','marker','.');
plot(wrap2_pert(1,1),lat2_pert(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap2_pert(end),lat2_pert(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - one day');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Perturbed','Start','End');
set(gca,'YDir','normal')

figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap3_pert,lat3_pert,'r','linestyle','none','marker','.');
plot(wrap3_pert(1,1),lat3_pert(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap3_pert(end),lat3_pert(end),'cs','MarkerSize',8,'LineWidth',2);
title('Ground Track - ten days');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Perturbed','Start','End');
set(gca,'YDir','normal')

k = 3;  m = 13;
[a_rep,T_rep] = RepeatGroundTrack(k,m,rad2deg(omegaE),mu_E);
tspan_rep = linspace(0, 2*T_rep, N);

% unperturbed repeating ground track
[alpha_rep,delta_rep,lon_rep,lat_rep] = ground_track(a_rep, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan_rep, omegaE, mu_E, t0);

lon_rep = rad2deg(lon_rep);
lat_rep = rad2deg(lat_rep);
wrap_rep = wrapTo180(lon_rep);

% Plot of unperturbed repeating ground tracks
figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap_rep,lat_rep,'g','linestyle','none','marker','.');
plot(wrap_rep(1,1),lat_rep(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap_rep(end),lat_rep(end),'cs','MarkerSize',8,'LineWidth',2);
title('Unperturbed Repeating Ground Track');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Unperturbed','Start','End');
set(gca,'YDir','normal')

% perturbed repeating ground track
[a_repJ2,T_repJ2] = RepeatGroundTrack_J2(k,m,omegaE,mu_E,e,i,a_rep);
tspan_repJ2 = linspace(0, 2*T_repJ2, N);

[alpha_rep_pert,delta_rep_pert,lon_rep_pert,lat_rep_pert] = ground_track_pert(a_repJ2, e, deg2rad(i), deg2rad(OM),...
    deg2rad(om), deg2rad(th) , thG_t0, tspan_repJ2, omegaE, mu_E, t0);

lon_rep_pert = rad2deg(lon_rep_pert);
lat_rep_pert = rad2deg(lat_rep_pert);
wrap_rep_pert = wrapTo180(lon_rep_pert);

% Plot of perturbed repeating ground tracks
figure()
img = imread('earthsurface.jpg');  % Load a sample image
image([-180 180],[90 -90],img);  % Plot the image
hold on
plot(wrap_rep_pert,lat_rep_pert,'r','linestyle','none','marker','.');
plot(wrap_rep_pert(1,1),lat_rep_pert(1,1),'ro','MarkerSize',8,'LineWidth',2);
plot(wrap_rep_pert(end),lat_rep_pert(end),'cs','MarkerSize',8,'LineWidth',2);
title('Perturbed Repeating Ground Track');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Perturbed','Start','End');
set(gca,'YDir','normal')


%% propagations

kep0 = [a, e, deg2rad(i), deg2rad(OM), deg2rad(om), deg2rad(th)];
[r0, v0] = kep2car([a,e,deg2rad(i),deg2rad(OM),deg2rad(om),deg2rad(th),mu_E]);
y0 = [r0; v0];
N = 10000;
initialdate = [2022,12,25,12,00,00];
initialmjd = date2mjd2000(initialdate);
t0 = initialmjd*(24*3600);
tspan = linspace(t0, t0+500*T, N);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Cartesian propagation:
[t_car, Y_car] = ode113(@(t,y) ode_2bp_perturbed(y,t),tspan,y0,options);
kep_car = zeros(N,6); 
for j=1:N                                                          
    kep_car(j,:) = car2kep(Y_car(j,1:3),Y_car(j,4:6),mu_E);
end
kep_car(:,4:6) = unwrap(kep_car(:,4:6));    

% Gauss propagation:
ode_Gauss = @(t,kep) GaussPert(kep,mu_E,t);
[t_Gauss,kep_gauss] = ode113(ode_Gauss,tspan,kep0,options);

% Errors between Gauss and Cartesian propagation:
err_a_plot = abs(kep_car(:,1)-kep_gauss(:,1))./abs(kep0(1));  % relative error of semi-major axis
err_e_plot = abs(kep_car(:,2)-kep_gauss(:,2));  % absolute error of eccentricity
err_i_plot = abs(kep_car(:,3)-kep_gauss(:,3))/(2*pi);  % relative error of inclination
err_OM_plot = abs(kep_car(:,4)-kep_gauss(:,4))/(2*pi);  % relative error of right ascension of ascending node
err_om_plot = abs(kep_car(:,5)-kep_gauss(:,5))/(2*pi);  % relative error of argoument of periapsis
err_th_plot = abs(kep_car(:,6)-kep_gauss(:,6))./abs(kep_gauss(:,6));  % relative error of true anomaly

figure()
subplot(2,3,1)
plot((t_car-t0)/T,kep_car(:,1),'LineWidth',2)
title('evolution of a _ Cartesian method [km]');
xlabel('n° orbit period');
ylabel('a_{car} [km]');
grid on

subplot(2,3,2)
plot((t_Gauss-t0)/T,kep_gauss(:,1),'LineWidth',2)
title('evolution of a _ Gauss method [km]');
xlabel('n° orbit period');
ylabel('a_{gauss} [km]');
grid on

subplot(2,3,3)
semilogy((t_Gauss-t0)/T,err_a_plot,'LineWidth',2)
title('rel. err. of a ');
xlabel('n° orbit period');
ylabel('|a_{car} - a_{gauss}|/|a_0| ');
grid on

subplot(2,3,4)
plot((t_car-t0)/T,kep_car(:,2),'LineWidth',2)
title('evolution of e _ Cartesian method [-]');
xlabel('n° orbit period');
ylabel('e_{car} [km]');
grid on

subplot(2,3,5)
plot((t_Gauss-t0)/T,kep_gauss(:,2),'LineWidth',2)
title('evolution of e _ Gauss method [-]');
xlabel('n° orbit period');
ylabel('e_{gauss} [-]');
grid on

subplot(2,3,6)
semilogy((t_Gauss-t0)/T,err_e_plot,'LineWidth',2)
title('abs. err. of e [-]');
xlabel('n° orbit period');
ylabel('|e_{car} - e_{gauss}|  [-]');
grid on


figure()
subplot(2,3,1)
plot((t_car-t0)/T,rad2deg(kep_car(:,3)),'LineWidth',2)
title('evolution of i _ Cartesian method [deg]');
xlabel('n° orbit period');
ylabel('i_{car} [deg]');
grid on

subplot(2,3,2)
plot((t_Gauss-t0)/T,rad2deg(kep_gauss(:,3)),'LineWidth',2)
title('evolution of i _ Gauss method [deg]');
xlabel('n° orbit period');
ylabel('i_{gauss} [deg]');
grid on

subplot(2,3,3)
semilogy((t_Gauss-t0)/T,rad2deg(err_i_plot),'LineWidth',2)
title('rel. err. of i');
xlabel('n° orbit period');
ylabel('|i_{car} - i_{gauss}|/360° ');
grid on

subplot(2,3,4)
plot((t_car-t0)/T,rad2deg(kep_car(:,4)),'LineWidth',2)
title('evolution of \Omega [deg] _ Cartesian method');
xlabel('n° orbit period');
ylabel('\Omega_{car} [deg]');
grid on

subplot(2,3,5)
plot((t_Gauss-t0)/T,rad2deg(kep_gauss(:,4)),'LineWidth',2)
title('evolution of \Omega [deg] _ Gauss method');
xlabel('n° orbit period');
ylabel('\Omega_{gauss} [deg]');
grid on

subplot(2,3,6)
semilogy((t_Gauss-t0)/T,rad2deg(err_OM_plot),'LineWidth',2)
title('rel. err. of \Omega ');
xlabel('n° orbit period');
ylabel('|\Omega_{car} - \Omega_{gauss}|/360°');
grid on
        
figure()
subplot(2,3,1)
plot((t_car-t0)/T,rad2deg(kep_car(:,5)),'LineWidth',2)
title('evolution of \omega _ Cartesian method [deg]');
xlabel('n° orbit period');
ylabel('\omega_{car} [deg]');
grid on

subplot(2,3,2)
plot((t_Gauss-t0)/T,rad2deg(kep_gauss(:,5)),'LineWidth',2)
title('evolution of \omega _ Gauss method [deg]');
xlabel('n° orbit period');
ylabel('\omega_{gauss} [deg]');
grid on

subplot(2,3,3)
semilogy((t_Gauss-t0)/T,rad2deg(err_om_plot),'LineWidth',2)
title('rel. err. of \omega');
xlabel('n° orbit period');
ylabel('|\omega_{car} - \omega_{gauss}|/360°');
grid on

subplot(2,3,4)
plot((t_car-t0)/T,rad2deg(kep_car(:,6)),'LineWidth',2)
title('evolution of \theta _ Cartesian method [deg]');
xlabel('n° orbit period');
ylabel('\theta_{car} [deg]');
grid on

subplot(2,3,5)
plot((t_Gauss-t0)/T,rad2deg(kep_gauss(:,6)),'LineWidth',2)
title('evolution of \theta _ Gauss method [deg]');
xlabel('n° orbit period');
ylabel('\theta_{gauss} [deg]');
grid on

subplot(2,3,6)
semilogy((t_Gauss-t0)/T,rad2deg(err_th_plot),'LineWidth',2)
title('rel. err. of \theta');
xlabel('n° orbit period');
ylabel('|\theta_{car} - \theta_{gauss}|/360°');
grid on

%% real data

% Initial keplerian elements of the satellite and other inputs:
kep0_sat = [4.046496437446377E+04, 2.917034269288031E-01, deg2rad(2.183757731769000E+01), deg2rad(2.616733747605937E+02), deg2rad(1.500679599533694E+02), deg2rad(1.372205091124491E+02)] ;

initial_date = [2019,01,01,00,00,00];           % initial date for timespan
final_date = [2020,01,01,00,00,00];             % final date for timespan
t_in = date2mjd2000(initial_date);              % initial date in MJD2000
t_fin = date2mjd2000(final_date);               % final date in MJD2000
t_in_sec = t_in*86400;                          % initial MJD2000 in seconds (for propagation purposes)
t_fin_sec = t_fin*86400;                        % initial MJD2000 in seconds (for propagation purposes)
Step = 86400;                                   % length in seconds of the time-step for propagation, chosen equal to the one of the ephemerides (1 day = 86400s)
t_span_sat = [t_in_sec:Step:t_fin_sec];         % defintion of the timespan

options = odeset('RelTol',1e-13,'AbsTol',1e-14);

% Propagation of the Ideal Prturbed Model 
GAUSS_PLANETARY_EQs = @(t,kep) ode_gauss_rsw_asgn (t,kep,mu_E,Keplerian_model_aj2_RSW(kep,mu_E,R_E,J2),TimeEph_model_a_moon_RSW(t,kep,mu_M,mu_E));
[Time_sat,kep_sat] = ode113(GAUSS_PLANETARY_EQs,t_span_sat,kep0_sat,options);

% Output of Real Ephemerides
namefileexcel= 'Real_ephemerides_ICRF-CUTforMATLAB.xlsx';       % excel file containg the ephemerides data             
[time_real_eph,kep_real_eph]= Real_ephemerides(namefileexcel);  % extracts the time and keplerian elements from the ephemerides file
T_sat= 2*pi*sqrt(kep0_sat(1)^3/mu_E);                           % orbital period of the satellite [s]

kep_real_eph(:,3:6) = deg2rad(kep_real_eph(:,3:6));             %converts keplerian elements in [rad]
kep_real_eph(:,4:6) = unwrap(kep_real_eph(:,4:6));              %unwraps keplerian elements obtained from the ephemerides
kep_sat(:,4:6) = unwrap(kep_sat(:,4:6));                        %unwraps keplerian elements obtained from the propagation

% Errors between Gauss and Cartesian propagation
err_a_plot1 = abs(kep_sat(:,1)-kep_real_eph(:,1))./abs(kep0_sat(1));         % relative error of semi-major axis
err_e_plot1 = abs(kep_sat(:,2)-kep_real_eph(:,2));                           % absolute error of eccentricity
err_i_plot1 = abs(kep_sat(:,3)-kep_real_eph(:,3))/(2*pi);                    % relative error of inclination
err_OM_plot1 = abs(kep_sat(:,4)-kep_real_eph(:,4))/(2*pi);                   % relative error of right ascension of ascending node
err_om_plot1 = abs(kep_sat(:,5)-kep_real_eph(:,5))/(2*pi);                   % relative error of argoument of periapsis
err_th_plot1 = abs(kep_sat(:,6)-kep_real_eph(:,6))./abs(kep_real_eph(:,6));  % relative error of true anomaly








































