function [R_pert,V_pert] = two_body_problem_fun_pert( r0, v0, mu, varargin )
%
% Two-Body Problem Solver for the perturbed problem.
%
% PROTOTYPE:
%   [R_pert,V_pert] = two_body_problem_fun( r0, v0, mu )
%
% DESCRIPTION
%   Solve the perturbed two body problem, integrating numerically the
%   equations, given the initial conditions, position and velocity.
%   Return the state vector of the body.
%
% INPUT:
%   r0 [3x1] Initial position of the body ( r0_x; r0_y; r0_z ) [km]
%   v0 [3x1] Initial velocity of the body ( v0_x; v0_y; v0_z ) [km/s]
%   varargin Variables to plot the orbit and it's main properties
%
% OUTPUT:
%   [R_pert,V_pert] State of the body [km, km/s]
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

set(0,'defaulttextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultlegendInterpreter','latex');

% Physical parameters
mu; % Gravitational parameter [km^3/s^2]
R_e = astroConstants(23); % Earth radius [km]

% Initial condition
r0; % Initial position [km]
v0; % Initial velocity [km/s]
y0 = [ r0; v0 ]; % Initial state vector

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu ); % Semi-major axis [km]
orbit_period = 2*pi*sqrt( a^3/mu ); % Orbital period [1/s]
year = 365*24*60*60; % year [s]
tspan = linspace( 0, 5*year, 500000 );

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
R = Y(:,1:3);
V = Y(:,4:6);

% Perform the integration for perturbated problem
[ T_pert, Y_pert ] = ode113( @(t,y) ode_2bp_perturbed( t, y), tspan, y0, options );
R_pert = Y_pert(:,1:3);
V_pert = Y_pert(:,4:6);


if nargin == 4

% Plot the orbit centered around the Earth
figure()
   scatter3( R_pert(:,1), R_pert(:,2), R_pert(:,3), 1, T_pert/orbit_period,'filled');  
   hold on ,
        [X,Y,Z] = sphere(50);
        X = -R_e*X;
        Y = -R_e*Y;
        Z = -R_e*Z;
        plot3( R(:,1), R(:,2), R(:,3), '.r','LineWidth',3);
        surf(gca,X,Y,Z,imread('earthsurface.jpg'),'LineStyle','none','FaceColor','texturemap');
        xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); 
        title('Perturbed Two - Body Problem orbit','FontName', 'Calibri','FontSize',14,'Color','k'); 
        legend('Perturbed orbit','Unperturbed orbit','Location','northeast'); grid on; axis equal;
        c = colorbar; c.Label.String = 'n° orbits';
    hold off

elseif nargin == 5

    % Constants of Integration
    h_pert = cross(R_pert,V_pert); % Angular Momentum [km^2/s]
    e_pert = cross(V_pert,h_pert)/mu - R_pert./vecnorm(R_pert,2,2); % Eccentricity vector [-]
    E_pert = (V_pert(:,1).^2 + V_pert(:,2).^2 + V_pert(:,3).^2)/2 - mu*(1./vecnorm(R_pert,2,2)); % Kinetic Energy [km^2/s^2]
    
    % Plot Angular Momentum in time
    figure()
    hold on
        title('Angular Momentum','FontWeight','bold','FontSize',14);
        plot(T_pert,vecnorm(h_pert,2,2),'.k','LineWidth',2,'DisplayName','$||h||$');
        plot(T_pert,h_pert(:,1),'.r','LineWidth',2,'DisplayName','$h_x$');
        plot(T_pert,h_pert(:,2),'.g','LineWidth',2,'DisplayName','$h_y$');
        plot(T_pert,h_pert(:,3),'.b','LineWidth',2,'DisplayName','$h_z$');
        ylabel('$||h|| ,\; h_x ,\; h_y ,\; h_z \; \left [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$');
        hl = legend('show');
        set(hl, 'Interpreter','latex')
        grid on;   
        hold off
    
    % Plot Eccentricity Vector in time
    figure()
    hold on
        title('Eccentricity Vector','FontWeight','bold','FontSize',14);
        plot(T_pert,vecnorm(e_pert,2,2),'.k','LineWidth',2,'DisplayName','$||e||$');
        plot(T_pert,e_pert(:,1),'.r','LineWidth',2,'DisplayName','$e_x$');
        plot(T_pert,e_pert(:,2),'.g','LineWidth',2,'DisplayName','$e_y$');
        plot(T_pert,e_pert(:,3),'.b','LineWidth',2,'DisplayName','$e_z$');
        ylabel('$||e||, \; e_x, \; e_y, \; e_z \; \left [ - \right ]$'); 
        xlabel('$t \left [ s \right ]$');
        hl = legend('show');
        set(hl, 'Interpreter','latex')
        grid on; 
    hold off
    
    % Plot dot product e-h in time
    figure()
    hold on
        title('$e \cdot h$','FontWeight','bold','FontSize',14);
        plot(T_pert,dot(e_pert,h_pert,2),'-b','LineWidth',2); 
        ylabel('$|e \cdot h| \; \left  [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$'); grid on;
    hold off
    
    % Plot Kinetic Energy in time
    figure()
    hold on
        title('Kinetic Energy','FontWeight','bold','FontSize',14);
        plot(T_pert,E_pert,'b','LineWidth',2); 
        ylabel('$E \; \left  [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$'); grid on;
    hold off

end

end