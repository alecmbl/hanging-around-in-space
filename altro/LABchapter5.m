clear
clc

%%

mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
R_e = astroConstants(23);
J2 = astroConstants(9);
N = 10000;

a = 7571;
e = 0.01;
i = deg2rad(87.9);
OM = deg2rad(180);
om = deg2rad(180);
th = deg2rad(0);
kep0 = [a, e, i, OM, om, th];

T = 2*pi*sqrt((a^3)/mu_E);
[r0, v0] = kep2car([a,e,i,OM,om,th,mu_E]);
y0 = [r0;v0];

T_orb = 2*pi*sqrt( a^3/mu_E ); % Orbital period [s]
tspan = linspace( 0, 100*T_orb, N);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Cartesian propagation:
[t_car, Y_car] = ode113(@(t,y) ode_2bp_J2(y),tspan,y0,options);
kep_car = zeros(N,6); 
for j=1:N                                                          
    kep_car(j,:) = car2kep(Y_car(j,1:3),Y_car(j,4:6),mu_E);
end
kep_car(:,4:6) = unwrap(kep_car(:,4:6)); 

% Gauss propagation:
ode_Gauss = @(t,kep) GaussPert_J2(t,kep,mu_E);
[t_Gauss,kep_gauss] = ode113(ode_Gauss,tspan,kep0,options);

% Errors between Gauss and Cartesian propagation:
err_a_plot = abs(kep_car(:,1)-kep_gauss(:,1))./abs(kep0(1));           % relative error of semi-major axis
err_e_plot = abs(kep_car(:,2)-kep_gauss(:,2));                          % absolute error of eccentricity
err_i_plot = abs(kep_car(:,3)-kep_gauss(:,3))/(2*pi);                   % relative error of inclination
err_OM_plot = abs(kep_car(:,4)-kep_gauss(:,4))/(2*pi);                  % relative error of right ascension of ascending node
err_om_plot = abs(kep_car(:,5)-kep_gauss(:,5))/(2*pi);                  % relative error of argoument of periapsis
err_th_plot = abs(kep_car(:,6)-kep_gauss(:,6))./abs(kep_gauss(:,6));    % relative error of true anomaly



%% plot

figure()
plot((t_car)/T,kep_car(:,1),'LineWidth',2)
title('evolution of a _ Cartesian method [km]');
xlabel('n° orbit period');
ylabel('a_{car} [km]');
grid on

figure()
plot((t_Gauss)/T,kep_gauss(:,1),'LineWidth',2)
title('evolution of a _ Gauss method [km]');
xlabel('n° orbit period');
ylabel('a_{gauss} [km]');
grid on

%%

figure()
plot((t_car)/T,kep_car(:,2),'LineWidth',2)
title('evolution of e _ Cartesian method [-]');
xlabel('n° orbit period');
ylabel('e_{car} [km]');
grid on

figure()
plot((t_Gauss)/T,kep_gauss(:,2),'LineWidth',2)
title('evolution of e _ Gauss method [-]');
xlabel('n° orbit period');
ylabel('e_{gauss} [-]');
grid on

%%

figure()
subplot(2,3,1)
plot((t_car)/T,kep_car(:,1),'LineWidth',2)
title('evolution of a _ Cartesian method [km]');
xlabel('n° orbit period');
ylabel('a_{car} [km]');
grid on

subplot(2,3,2)
plot((t_Gauss)/T,kep_gauss(:,1),'LineWidth',2)
title('evolution of a _ Gauss method [km]');
xlabel('n° orbit period');
ylabel('a_{gauss} [km]');
grid on

subplot(2,3,3)
semilogy((t_Gauss)/T,err_a_plot,'LineWidth',2)
title('rel. err. of a ');
xlabel('n° orbit period');
ylabel('|a_{car} - a_{gauss}|/|a_0| ');
grid on

subplot(2,3,4)
plot((t_car)/T,kep_car(:,2),'LineWidth',2)
title('evolution of e _ Cartesian method [-]');
xlabel('n° orbit period');
ylabel('e_{car} [km]');
grid on

subplot(2,3,5)
plot((t_Gauss)/T,kep_gauss(:,2),'LineWidth',2)
title('evolution of e _ Gauss method [-]');
xlabel('n° orbit period');
ylabel('e_{gauss} [-]');
grid on

subplot(2,3,6)
semilogy((t_Gauss)/T,err_e_plot,'LineWidth',2)
title('abs. err. of e [-]');
xlabel('n° orbit period');
ylabel('|e_{car} - e_{gauss}|  [-]');
grid on





















