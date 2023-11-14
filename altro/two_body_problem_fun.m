function [R,V] = two_body_problem_fun( r0, v0, mu, varargin )
%
% Two-Body Problem Solver
%
% PROTOTYPE:
%   [R,V] = two_body_problem_fun( r0, v0, mu )
%
% DESCRIPTION
%   Solve the unperturbed two body problem, integrating numerically the
%   equations, given the initial conditions, position and velocity.
%   Return the state vector of the body fo one orbital period.
%
% INPUT:
%   r0[3x1] Initial position of the body ( r0_x; r0_y; r0_z ) [km]
%   v0[3x1] Initial velocity of the body ( v0_x; v0_y; v0_z ) [km/s]
%   varargin Variables to plot the orbit and it's main properties
%
% OUTPUT:
%   [R,V] State of the body for one orbit period [km, km/s]
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

% Initial condition
r0; % Initial position [km]
v0; % Initial velocity [km/s]
y0 = [ r0; v0 ]; % Initial state vector

% Set time span
a = 1/( 2/norm(r0) - dot(v0,v0)/mu ); % Semi-major axis [km]
orbit_period = 2*pi*sqrt( a^3/mu ); % Orbital period [1/s]
tspan = linspace( 0, orbit_period, 1000 ); % Time span

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );
R = Y(:,1:3);
V = Y(:,4:6);

% Plot the orbit centered about the earth
if nargin == 4

    figure()
        plot3( R(:,1), R(:,2), R(:,3), '.r' );
        xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]'); grid on; axis equal;
        title('Two-body problem orbit','FontWeight','bold','FontSize',14,'Color','k');
        hold on
            h1= gca;
            earth_sphere(h1,'km');
        hold off   

elseif nargin == 5

    % Constants of Integration
    h = cross(R,V); % Angular Momentum [km^2/s]
    e = cross(V,h)/mu - R./vecnorm(R,2,2); % Eccentricity vector [-]
    E = (V(:,1).^2 + V(:,2).^2 + V(:,3).^2)/2 - mu*(1./vecnorm(R,2,2)); % Kinetic Energy [km^2/s^2]
    
    
    %Plot Angular Momentum in time
    figure()
    hold on
        title('Angular Momentum','FontWeight','bold','FontSize',14);
        plot(T,vecnorm(h,2,2),'--k','LineWidth',3,'DisplayName','$||h||$');
        plot(T,h(:,1),'--r','LineWidth',3,'DisplayName','$h_x$');
        plot(T,h(:,2),'--g','LineWidth',3,'DisplayName','$h_y$');
        plot(T,h(:,3),'--b','LineWidth',3,'DisplayName','$h_z$');
        ylabel('$||h|| ,\; h_x ,\; h_y ,\; h_z \; \left [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$');
        hl = legend('show');
        set(hl, 'Interpreter','latex')
        grid on; ylim([-11e4 11e4]);   
        hold off
    
    %Plot Eccentricity Vector in time
    figure()
    hold on
        title('Eccentricity Vector','FontWeight','bold','FontSize',14);
        plot(T,vecnorm(e,2,2),'--k','LineWidth',3,'DisplayName','$||e||$');
        plot(T,e(:,1),'--r','LineWidth',1.75,'DisplayName','$e_x$');
        plot(T,e(:,2),'--g','LineWidth',3,'DisplayName','$e_y$');
        plot(T,e(:,3),'--b','LineWidth',1.75,'DisplayName','$e_z$');
        ylabel('$||e||, \; e_x, \; e_y, \; e_z \; \left [ - \right ]$'); 
        xlabel('$t \left [ s \right ]$');
        hl = legend('show');
        set(hl, 'Interpreter','latex')
        grid on; ylim([-0.5e-4 2.5e-4]); 
    hold off
    
    %Plot dot product e-h in time
    figure()
    hold on
        title('$e \cdot h$','FontWeight','bold','FontSize',14);
        %e_dot_p = e(:,1).*h(:,1)+e(:,1).*h(:,1)+e(:,1).*h(:,1);
        plot(T,dot(e,h,2),'-b','LineWidth',2); 
        ylabel('$|e \cdot h| \; \left  [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$'); grid on;
    hold off
    
    %Plot Kinetic Energy in time
    figure()
    hold on
        title('Kinetic Energy','FontWeight','bold','FontSize',14);
        plot(T,E,'b','LineWidth',2); 
        ylabel('$E \; \left  [ km^2/s \right ]$'); 
        xlabel('$t \left [ s \right ]$'); grid on;
hold off
end

end

