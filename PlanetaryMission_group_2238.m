
% This is the code for the Assignment 2. 
% The script gives all the results required:
%   - Orbit Plot (unperturbed and perturbed)
%   - Ground Tracks for all assigned conditions
%   - Gauss propagation of Keplerian elements
%   - Keplerian Element Errors between Cartesian and Gauss propagation
%   - Comparison with real data from Ephemerides

%%

clc;
clear all;
close all;

%% Path

addpath('functions')
addpath('Default_function')
addpath('time')
set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultlegendInterpreter','latex');

%% Initial conditions

mu_E = astroConstants(13);             % Earth's gravitational parameter [km^3/s^2]
R_e = astroConstants(23);              % Earth's radius [km]
J2 = astroConstants(9);                % J2 constant

thG_t0 = 0;                            % theta of Greenwich Meridian [rad]
omegaE = ((15.04*pi)/180)/3600;        % Earth Angular Velocity [rad/s]

t0 = 0;                                % Initial time [s]
th0 = 0;                               % Initial True Anomaly [rad]
i_rad = (46.2887*pi)/180;              % Inclination [rad]
a = 109394;                            % Semi-Major axis [km]
e = 0.3100;                            % Eccentricity [-]
OM = 0.0087;                           % Right Ascention of ascending node [rad ]
om = 0;                                % Argument of the Pericenter [rad]
T = 2*pi*sqrt( a^3/mu_E );             % Orbital period [s]
kep_0 = [a,e,i_rad,OM,om,th0,mu_E];    % Initial keplerian elements
[r0, v0] = kep2car(kep_0);             
y0 = [ r0; v0 ];                       % Initial cartesian state

time_string = ['1 Orbit';'1 Day  ';'10 Days'];


%% Nominal Orbit Plot

[R,V] = two_body_problem_fun( r0, v0, mu_E, 1 );


%% Perturbed Orbit Plot
[R_pert,V_pert] = two_body_problem_fun_pert( r0, v0, mu_E, 1 );

% cartesian state and keplerian elements after perturbed propagation
Final_state = [R_pert(end,1); R_pert(end,2);R_pert(end,3)];
Final_velocity = [V_pert(end,1); V_pert(end,2);V_pert(end,3)];
Final_kep = car2kep(Final_state,Final_velocity,mu_E); 
T_final = 2*pi*sqrt( Final_kep(1)^3/mu_E);

%% Ground Track (no perturbations)

% Set time span
t_vect = [linspace( 0, T, 10000 );linspace( 0, 24*3600, 10000 );linspace( 0, 10*24*3600, 10000 )];
% t_vect(1) --> 1 orbit, t_vect(2) --> 1 day, t_vect(3) --> 10 days

for count = 1:3

    [alpha,delta,lon,lat] = ground_track(a, e, i_rad, OM, om, th0 , thG_t0, t_vect(count,:), omegaE, mu_E, t0);

    % conversion from radiants to degrees
    lon_deg = rad2deg(lon); 
    lat_deg = rad2deg(lat);
    wrap = wrapTo180(lon_deg);
    
    % plot
    figure()
    img = imread('earthsurface.jpg');  % Load a sample image
    image([-180 180],[90 -90],img);  % Plot the image
    hold on
    plot(wrap,lat_deg,'g','LineStyle', 'none', 'Marker', '.')
    hold on
    plot(wrap(1,1),lat_deg(1,1),'o',wrap(end),lat_deg(end),'o','LineWidth',2)
    xlim([-180 180]); ylim([-90 90]); grid on;
    xlabel('$Longitude \left [ deg \right ]$'); ylabel('$Latitude \left [ deg \right ]$');
    xticks(-180:30:180);
    yticks(-90:30:90);
    title(sprintf('Ground Track for %s',string(time_string(count,:))),'FontWeight','bold','FontSize',14,'Color','k');
    legend('Ground Track', 'Start', 'End');
    set(gca,'YDir','normal')

end

%% Repeating Ground Track (no perturbations)

k = 3; % number of satellite revolutions
m = 13; % number of earth revolutions

[a_rep,T_rep] = RepeatGroundTrack(k,m,rad2deg(omegaE),mu_E);

% set time span
t_rep = [linspace( 0, T_rep, 10000 );linspace( 0, 24*3600, 10000 );linspace( 0, 10*24*3600, 10000 )];
% t_rep(1) --> 1 orbit, t_rep(2) --> 1 day, t_rep(3) --> 10 days

count = 0;

for count = 1:3

    [alpha_rep,delta_rep,lon_rep,lat_rep] = ground_track(a_rep, e, i_rad, OM, om, th0 , thG_t0, t_rep(count,:), omegaE, mu_E, t0);

    % conversion from radiants to degrees
    lon_deg_2 = rad2deg(lon_rep);
    lat_deg_2 = rad2deg(lat_rep);
    wrap_rep = wrapTo180(lon_deg_2);

    % plot
    figure()
    img = imread('earthsurface.jpg');  % Load a sample image
    image([-180 180],[90 -90],img);  % Plot the image
    hold on
    plot(wrap_rep,lat_deg_2,'r','linestyle','none','marker','.')
    hold on
    plot(wrap_rep(1,1),lat_deg_2(1,1),'o',wrap_rep(end),lat_deg_2(end),'o','LineWidth',2)
    xlim([-180 180]); ylim([-90 90]); grid on;
    xlabel('$Longitude \left [ deg \right ]$'); ylabel('$Latitude \left [ deg \right ]$'); 
    xticks(-180:30:180);
    yticks(-90:30:90);
    title(sprintf('Repeating Ground Track for %s',string(time_string(count,:))),'FontWeight','bold','FontSize',14,'Color','k');
    legend('Repeating Ground Track', 'Start', 'End');
    set(gca,'YDir','normal')

end

%% Ground Track (perturbations)

% Set time span
t_vect = [linspace( 0, T, 10000 );linspace( 0, 24*3600, 10000 );linspace( 0, 10*24*3600, 10000 )];
% t_vect(1) --> 1 orbit, t_vect(2) --> 1 day, t_vect(3) --> 10 days

for count = 1:3

    [alpha,delta,lon,lat] = ground_track_perturbed(a, e, i_rad, OM, om, th0 , thG_t0, t_vect(count,:), omegaE, mu_E, t0);

    % conversion from radiants to degrees
    lon_deg = rad2deg(lon); 
    lat_deg = rad2deg(lat);
    wrap = wrapTo180(lon_deg);
    
    % plot
    figure()
    img = imread('earthsurface.jpg');  % Load a sample image
    image([-180 180],[90 -90],img);  % Plot the image
    hold on
    plot(wrap,lat_deg,'g','LineStyle', 'none', 'Marker', '.')
    hold on
    plot(wrap(1,1),lat_deg(1,1),'o',wrap(end),lat_deg(end),'o','LineWidth',2)
    xlim([-180 180]); ylim([-90 90]); grid on;
    xlabel('$Longitude \left [ deg \right ]$'); ylabel('$Latitude \left [ deg \right ]$'); 
    xticks(-180:30:180);
    yticks(-90:30:90);
    title(sprintf('Perturbed Ground Track for %s',string(time_string(count,:))),'FontWeight','bold','FontSize',14,'Color','k');
    legend('Ground Track', 'Start', 'End');
    set(gca,'YDir','normal')

end

%% Repeating Ground Track (perturbations)

k = 3; % number of satellite revolutions
m = 13; % number of earth revolutions

[a_repJ2,T_repJ2] = RepeatGroundTrack_J2(k,m,omegaE,mu_E,e,i_rad,a_rep);

% set time span
t_rep = [linspace( 0, T_repJ2, 10000 );linspace( 0, 24*3600, 10000 );linspace( 0, 10*24*3600, 10000 )];
% t_rep(1) --> 1 orbit, t_rep(2) --> 1 day, t_rep(3) --> 10 days

count = 0;

for count = 1:3

    [alpha_rep,delta_rep,lon_rep,lat_rep] = ground_track_perturbed(a_repJ2, e, i_rad, OM, om, th0 , thG_t0, t_rep(count,:), omegaE, mu_E, t0);

    % conversion from radiants to degrees
    lon_deg_2 = rad2deg(lon_rep);
    lat_deg_2 = rad2deg(lat_rep);
    wrap_rep = wrapTo180(lon_deg_2);

    % plot
    figure()
    img = imread('earthsurface.jpg');  % Load a sample image
    image([-180 180],[90 -90],img);  % Plot the image
    hold on
    plot(wrap_rep,lat_deg_2,'r','linestyle','none','marker','.')
    hold on
    plot(wrap_rep(1,1),lat_deg_2(1,1),'o',wrap_rep(end),lat_deg_2(end),'o','LineWidth',2)
    xlim([-180 180]); ylim([-90 90]); grid on;
    xlabel('$Longitude \left [ deg \right ]$'); ylabel('$Latitude \left [ deg \right ]$'); 
    xticks(-180:30:180);
    yticks(-90:30:90);
    title(sprintf('Perturbed Repeating Ground Track for %s',string(time_string(count,:))),'FontWeight','bold','FontSize',14,'Color','k');
    legend('Repeating Ground Track', 'Start', 'End');
    set(gca,'YDir','normal')
end

%% Cartesian - Gauss Propagations 

kep0 = [a, e, i_rad, OM, om, th0];
[r0, v0] = kep2car([a,e,i_rad,OM,om,th0,mu_E]);
y0 = [r0; v0];
N = 1000;
initialdate = [2022,12,25,12,00,00];
initialmjd = date2mjd2000(initialdate);
t0 = initialmjd*(24*3600);
periods = 3; % n° of orbit perturbation period
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Cartesian propagation:
[t_car, Y_car] = ode113(@(t,y) ode_2bp_perturbed(t,y),linspace(t0, t0 + periods*30*24*3600, N),y0,options);
kep_car = zeros(N,6); 
for j=1:N                                                          
    kep_car(j,:) = car2kep(Y_car(j,1:3),Y_car(j,4:6),mu_E);
end
kep_car(:,4:6) = unwrap(kep_car(:,4:6));  
t_plot_car = (t_car-t0)/(24*3600);

% Gauss propagation:
ode_Gauss = @(t,kep) GaussPert(kep,mu_E,t);
[t_Gauss,kep_gauss] = ode113(ode_Gauss,linspace(t0, t0 + periods*30*24*3600, N),kep0,options);
t_plot_Gauss = (t_Gauss-t0)/(24*3600);

% Errors between Gauss and Cartesian propagation:
err_a_plot = abs(kep_car(:,1)-kep_gauss(:,1))./abs(kep0(1));            % relative error of semi-major axis
err_e_plot = abs(kep_car(:,2)-kep_gauss(:,2));                          % absolute error of eccentricity
err_i_plot = abs(kep_car(:,3)-kep_gauss(:,3))/(2*pi);                   % relative error of inclination
err_OM_plot = abs(kep_car(:,4)-kep_gauss(:,4))/(2*pi);                  % relative error of right ascension of ascending node
err_om_plot = abs(kep_car(:,5)-kep_gauss(:,5))/(2*pi);                  % relative error of argoument of periapsis
err_th_plot = abs(kep_car(:,6)-kep_gauss(:,6))./abs(kep_gauss(:,6));    % relative error of true anomaly

%% Propagation Plots New 

figure()
plot(t_plot_car,kep_car(:,1),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('Evolution of a (Cartesian) [km]');
xlabel('$n ^\circ \;days$');
ylabel('$a_{car}$ [km]');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,kep_gauss(:,1),'LineWidth',2,'Color','r')
title('Evolution of a (Gauss) [km]');
xlabel('$n ^\circ \;days$');
ylabel('$a_{gauss}$ [km]');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,err_a_plot,'LineWidth',2,'Color','k')
title('Rel. Err. of a ');
xlabel('$n ^\circ \;days$');
ylabel('$|a_{car} - a_{gauss}|/|a_0|$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_car,kep_car(:,2),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('Evolution of e (Cartesian) [-]');
xlabel('$n ^\circ \;days$');
ylabel('$e_{car}$ [-]');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,kep_gauss(:,2),'LineWidth',2,'Color','r')
title('Evolution of e (Gauss) [-]');
xlabel('$n ^\circ \;days$');
ylabel('$e_{gauss}$ [-]');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,err_e_plot,'LineWidth',2,'Color','k')
title('Abs. Err. of e [-]');
xlabel('$n ^\circ \;days$');
ylabel('$|e_{car} - e_{gauss}| [-]$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_car,rad2deg(kep_car(:,3)),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('Evolution of i (Cartesian) [deg]');
xlabel('$n ^\circ \;days$');
ylabel('$i_{car}$ [deg]');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,rad2deg(kep_gauss(:,3)),'LineWidth',2,'Color','r')
title('Evolution of i (Gauss) [deg]');
xlabel('$n ^\circ \;days$');
ylabel('$i_{gauss}$ [deg]');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,rad2deg(err_i_plot),'LineWidth',2,'Color','k')
title('Rel. Err. of i');
xlabel('$n ^\circ \;days$');
ylabel('$|i_{car} - i_{gauss}|/360^\circ$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_car,rad2deg(kep_car(:,4)),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('$Evolution \;of \; \Omega\; (Cartesian) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\Omega_{car} [deg]$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,rad2deg(kep_gauss(:,4)),'LineWidth',2,'Color','r')
title('$Evolution\; of\; \Omega\;(Gauss) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\Omega_{gauss} [deg]$');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,rad2deg(err_OM_plot),'LineWidth',2,'Color','k')
title('$Rel.\; Err.\; of \;\Omega$');
xlabel('$n ^\circ \;days$');
ylabel('$|\Omega_{car} - \Omega_{gauss}|/360^\circ$');
xlim([0 periods*30]);
grid on
        
figure()
plot(t_plot_car,rad2deg(kep_car(:,5)),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('$Evolution\; of\; \omega\; (Cartesian) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\omega_{car} [deg]$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,rad2deg(kep_gauss(:,5)),'LineWidth',2,'Color','r')
title('$Evolution\; of\; \omega\; (Gauss) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\omega_{gauss} [deg]$');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,rad2deg(err_om_plot),'LineWidth',2,'Color','k')
title('$Rel.\; Err.\; of\; \omega$');
xlabel('$n ^\circ \;days$');
ylabel('$|\omega_{car} - \omega_{gauss}|/360^\circ$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_car,rad2deg(kep_car(:,6)),'LineWidth',2,'Color',[0 0.4470 0.7410])
title('$Evolution\; of\; \theta\;(Cartesian) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\theta_{car} [deg]$');
xlim([0 periods*30]);
grid on

figure()
plot(t_plot_Gauss,rad2deg(kep_gauss(:,6)),'LineWidth',2,'Color','r')
title('$Evolution\; of\; \theta\; (Gauss) [deg]$');
xlabel('$n ^\circ \;days$');
ylabel('$\theta_{gauss} [deg]$');
xlim([0 periods*30]);
grid on

figure()
semilogy(t_plot_Gauss,rad2deg(err_th_plot),'LineWidth',2,'Color','k')
title('$Rel.\; Err.\; of\; \theta$');
xlabel('$n ^\circ \;days$');
ylabel('$|\theta_{car} - \theta_{gauss}|/360^\circ$');
xlim([0 periods*30]);
grid on;
 
%% Filtering

% a
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,1)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,1)),50,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,1)),10000,'Endpoints','shrink'),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; a\; filtered$','FontWeight','bold','FontSize',12); 
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('a [km]');
hold off

% e
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,2)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,2)),70,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,2)),10000),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; e\; filtered$','FontWeight','bold','FontSize',12);  
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('e [-]');
hold off
 
% i
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,3)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,3)),70,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    coefficents = polyfit([t_plot_Gauss(1) t_plot_Gauss(end)],[rad2deg(kep_gauss(1,3)) rad2deg(kep_gauss(end,3))],1);
    e_lin = @ (x) coefficents(1).*x + coefficents(2);
    plot(t_plot_Gauss, e_lin(t_plot_Gauss),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; i\; filtered$','FontWeight','bold','FontSize',12); 
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('i [deg]');
hold off

% Omega
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,4)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,4)),50,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    coefficents = polyfit([t_plot_Gauss(1) t_plot_Gauss(end)],[rad2deg(kep_gauss(1,4)) rad2deg(kep_gauss(end,4))],1);
    Omega_lin = @ (x) coefficents(1).*x + coefficents(2);
    plot(t_plot_Gauss, Omega_lin(t_plot_Gauss),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; \Omega\; filtered$','FontWeight','bold','FontSize',12); 
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('$\Omega$ [deg]');
hold off

% omega
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,5)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,5)),70,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    coefficents = polyfit([t_plot_Gauss(1) t_plot_Gauss(end)],[rad2deg(kep_gauss(1,5)) rad2deg(kep_gauss(end,5))],1);
    omega_lin = @ (x) coefficents(1).*x + coefficents(2);
    plot(t_plot_Gauss, omega_lin(t_plot_Gauss),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; \omega\; filtered$','FontWeight','bold','FontSize',12); 
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('$\omega$ [deg]');
hold off

% theta
figure()
hold on
    plot(t_plot_Gauss,rad2deg(kep_car(:,6)),'DisplayName','unfiltered','Color',[0 0.4470 0.7410],'LineWidth',1.2);
    plot(t_plot_Gauss,movmean(rad2deg(kep_gauss(:,6)),10,'Endpoints','shrink'),'DisplayName','filtered','Color','r','LineWidth',1.2);
    coefficents = polyfit([t_plot_Gauss(1) t_plot_Gauss(end)],[rad2deg(kep_gauss(1,6)) rad2deg(kep_gauss(end,6))],1);
    theta_lin = @ (x) coefficents(1).*x + coefficents(2);
    plot(t_plot_Gauss, theta_lin(t_plot_Gauss),'DisplayName','secular','Color',[0.9290 0.6940 0.1250],'LineWidth',1.2);
    title('$Evolution\; of\; \theta\; filtered$','FontWeight','bold','FontSize',12); 
    subtitle('(Gauss method)','FontSize',10);
    grid on; xlim([0 periods*30]); legend(); xlabel('$n ^\circ \;days$'); ylabel('$\theta$ [deg]');
hold off

%% Comparison with real data from Ephemerides

% satellite GEOTAIL 22049 (Norad-ID)

% load real data from ephemerides
[kep_real, EpochY, EpochD] = readTLE('spacetrack_results.txt');

% from initial real kep elements to initial real car elements
kep_real_0 = kep_real(1,:);
[r0, v0] = kep2car([kep_real_0,mu_E]);
y0 = [r0; v0];

% changing time reference
time = [];
for j=1:1:length(EpochY)
      year = 2000 + EpochY(j);
      date = [year 01 01 00 00 00];
      time = [time; date2mjd2000(date) + EpochD(j)];
end

% set time span
t0 = time(1)*24*3600;
tf = time(end)*24*3600;

% cartesian propagation
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t_real_prop, Y_real_prop] = ode113(@(t,y) ode_2bp_perturbed(t,y),linspace(t0,tf,size(kep_real,1)),y0,options);
t_real_prop = t_real_prop./(24*3600) - time(1);

% from car to kep elements propagated
kep_real_prop = zeros(size(Y_real_prop,1),6);
for j = 1:size(Y_real_prop,1)                                                        
    kep_real_prop(j,:) = car2kep(Y_real_prop(j,1:3),Y_real_prop(j,4:6),mu_E);
end

%% Plot real data orbit from Ephemerides vs. progated orbit

figure;
plot(t_real_prop,kep_real(:,1));
hold on 
    plot(t_real_prop,kep_real_prop(:,1));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$a \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southwest');
    xlabel('$n ^\circ \;days$'); 
    ylabel('a [km]');
hold off

figure;
plot(t_real_prop,kep_real(:,2));
hold on 
    plot(t_real_prop,kep_real_prop(:,2));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$e \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southwest');
    xlabel('$n ^\circ \;days$');
    ylabel('e [-]');
hold off

figure;
plot(t_real_prop,unwrap(rad2deg(kep_real(:,3))));
hold on 
    plot(t_real_prop,unwrap(rad2deg(kep_real_prop(:,3))));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$i \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southwest');
    xlabel('$n ^\circ \;days$');
    ylabel('i [deg]');
hold off

figure;
plot(t_real_prop,unwrap(rad2deg(kep_real(:,4))));
hold on 
    plot(t_real_prop,unwrap(rad2deg(kep_real_prop(:,4))));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$\Omega \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southwest');
    xlabel('$n ^\circ \;days$');
    ylabel('$\Omega$ [deg]');
hold off

figure;
plot(t_real_prop,unwrap(rad2deg(kep_real(:,5))));
hold on 
    plot(t_real_prop,unwrap(rad2deg(kep_real_prop(:,5))));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$\omega \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southeast');
    xlabel('$n ^\circ \;days$');
    ylabel('$\omega$ [deg]');
hold off

figure;
plot(t_real_prop,unwrap(rad2deg(kep_real(:,6))));
hold on 
    plot(t_real_prop,unwrap(rad2deg(kep_real_prop(:,6))));
    xline([1090 2850],'-',{'Manoeuvre Begin'});
    xline([1500 3000],'-',{'Manoeuvre End'});
    title('$\theta \; comparison \; Ephemerides - Propagated \; data$','FontSize',12);
    xlim([0 3500]); grid on; 
    legend('Ephemerides','Propagated','Location','southwest');
    xlabel('$n ^\circ \;days$');
    ylabel('$\theta$ [deg]');
hold off
