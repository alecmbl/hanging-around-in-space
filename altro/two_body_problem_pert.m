close all
clear all
clc

%% Perturbed Two Body Problem

% Physical parameters
mu_E = astroConstants(13); % Earth's gravitational parameter [km^3/s^2]
J2_E = astroConstants(9); % Perturbation Constant [-]
R_E = astroConstants(23); % Earth Radius [km]

% Initial condition
t0 = 0;
th0 = 0;
i_rad = (46.2887*pi)/180; % [rad]
a = 109394*1e3; %[km]
e = 0.3100;
OM = 0;
om = 0;
kep = [109394,0.3100,i_rad,0,0,0,mu_E];
[r, v] = kep2car(kep);
y0 = [ r; v ];

% Set time span
orbit_period = 2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
month = 60 * 60 * 24 * 30;
year = 365*24*60*60; % year [s]
tspan = linspace( 0, 2000*month, 10000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration for perturbater problem
[ T_pert, Y_pert ] = ode113( @(t,y) ode_2bp_perturbed(t, y), tspan, y0, options );
R_pert = Y_pert(:,1:3);
V_pert = Y_pert(:,4:6);
% Perform the integration for base problem
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
R = Y(:,1:3);
V = Y(:,4:6);


%Plot the orbits
figure()
    scatter3( R_pert(:,1), R_pert(:,2), R_pert(:,3), 1, T_pert,'DisplayName','perturbed problem (J2)'); 
    hold on 
        plot3( R(:,1), R(:,2), R(:,3), '-r','LineWidth',1,'DisplayName','unperturbed problem (2BP)' );
        xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');  
        legend('show','Location','northeast');grid on; axis equal;colorbar;
        title('Perturbed (J2) Two-Body Problem orbit','FontName', 'Calibri','FontSize',14,'Color','k');
    hold off
%%


% Quantities
h_pert = cross(R_pert,V_pert); %Angular Momentum
e_pert = cross(V_pert,h_pert)/mu_E - R_pert./vecnorm(R_pert,2,2); %Eccentricity vector
E_pert = (V_pert(:,1).^2 + V_pert(:,2).^2 + V_pert(:,3).^2)/2 - mu_E*(1./vecnorm(R_pert,2,2)); %Kinetic Energy
h = cross(R,V); %Angular Momentum
e = cross(V,h)/mu_E - R./vecnorm(R,2,2); %Eccentricity vector
E = (V(:,1).^2 + V(:,2).^2 + V(:,3).^2)/2 - mu_E*(1./vecnorm(R,2,2)); %Kinetic Energy


T_pert = T_pert./(60*60*24*30);
T = T./(60*60*24*30);

%Plot Angular Momentum
figure()
    hold on
        title('Angular Momentum');
        plot(T_pert,vecnorm(h_pert,2,2),'--k','LineWidth',1.75,'DisplayName','||h||, J2');
        plot(T_pert,h_pert(:,1),'--r','LineWidth',1.75,'DisplayName','h_x, J2');
        plot(T_pert,h_pert(:,2),'--g','LineWidth',1.75,'DisplayName','h_y, J2');
        plot(T_pert,h_pert(:,3),'--b','LineWidth',1.75,'DisplayName','h_z, J2');
        plot(T,vecnorm(h,2,2),'--k','LineWidth',2.5,'DisplayName','||h||, 2BP');
        plot(T,h(:,1),'--r','LineWidth',2.5,'DisplayName','h_x, 2BP');
        plot(T,h(:,2),'--g','LineWidth',2.5,'DisplayName','h_z, 2BP');
        ylabel('[km^2/s]'); xlabel('t [s]');
        legend ('show'); grid on; ylim([-11e4 11e4]);   
    hold off
%%
%Plot Eccentricity Vector
figure()
    hold on
        title('Eccentricity vector');
        plot(T_pert,vecnorm(e_pert,2,2),'k','LineWidth',1.75,'DisplayName','||e||, J2');
        plot(T_pert,e_pert(:,1),'r','LineWidth',1.75,'DisplayName','e_x, J2');
        plot(T_pert,e_pert(:,2),'g','LineWidth',1.75,'DisplayName','e_y, J2');
        plot(T_pert,e_pert(:,3),'b','LineWidth',1.75,'DisplayName','e_z, J2');
        plot(T,vecnorm(e,2,2),'k','LineWidth',2.5,'DisplayName','\lVert e \rVert 2BP');
        plot(T,e(:,1),'r','LineWidth',2.5,'DisplayName','e_x, 2BP');
        plot(T,e(:,2),'g','LineWidth',2.5,'DisplayName','e_y, 2BP');
        plot(T,e(:,3),'b','LineWidth',2.5,'DisplayName','e_z, 2BP');
        ylabel('[-]'); xlabel('t [s]');
        legend ('show'); grid on; 
    hold off
%%
%Plot dot product e-h
figure()
    hold on
        title('e dot h');
        plot(T,dot(e,h,2),'-b','DisplayName','2BP','LineWidth',3);
        plot(T_pert,dot(e_pert,h_pert,2),'-r','DisplayName','J2-Perturbed','LineWidth',2); 
        ylabel('[km^2/s]'); xlabel('t [s]'); legend ('show'); grid on;
    hold off

%Plot Kinetic Energy
figure()
    hold on
        title('Kinetic energy');
        plot(T,E,'-b','DisplayName','2BP','LineWidth',3);
        plot(T_pert,E_pert,'--r','DisplayName','J2-Perturbed','LineWidth',2);
        ylabel('[km^2/s^2]'); xlabel('t [s]'); legend ('show'); grid on;
    hold off


