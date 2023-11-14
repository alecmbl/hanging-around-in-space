function plot_powerGA(dates, DepartureID, FlybyID, ArrivalID)
%
% Plot of the hyperbola of a flyby in the plante's perifocal frame
% 
% PROTOTYPE:
%   plot_powerGA(dates, DepartureID, FlybyID, ArrivalID)
%
% DESCRIPTION:
%   Plots the incoming and outgoing hyperbolas of a fly-by manouver
%   together with the powered fraction of delta V. The vector is scaled by
%   a factor of 1e8 to make it visible.
% 
% INPUT:
%   dates             [3x1]   Departure-FlyBy-Arrival dates      [datetime]
%   DepartureID       [1x1]   Departure planet ID                [-]
%   FlybyID           [1x1]   Fly-by planet ID                   [-]
%   ArrivalID         [1x1]   Arrival planet ID                  [-]
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

%% Solver options

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

%% Date format

date_vector = strcmp(class(dates(1)), 'double');
date_time = strcmp(class(dates(2)), 'datetime');

if date_vector
    Dep_date = datetime(dates(1));
    Arr_date = datetime(dates(3));
    Fb_date = datetime(dates(2));

    [Y, M, D] = ymd(datetime(dates(2)));
    [H, m, s] = hms(datetime(dates(2)));
    GA_date = [Y M D H m s];

    t_GA = date2mjd2000(GA_date);     % [days]

elseif date_time
    Dep_date = dates(1);
    Arr_date = dates(3);
    Fb_date = dates(2);

    t_GA = date2mjd2000(datevec(dates(2)));  % [days]

end

%% Constants

AU = astroConstants(2);
ksun = astroConstants(4);

%% Initial and final velocities

% Transfer arc 1
[~, ~, ~, ~, V_helio_inFb] = lambert_transfer(Dep_date, Fb_date,...
    seconds(Fb_date - Dep_date), DepartureID, FlybyID, ksun, 0);

% Transfer arc 2
[~, ~, ~, V_helio_outFb, ~] = lambert_transfer(Fb_date, Arr_date,...
    seconds(Arr_date - Fb_date), FlybyID, ArrivalID, ksun, 0);

%% Planet

[kep, ~] = uplanet(t_GA, FlybyID);
mu = astroConstants(FlybyID + 10);

[~, v_p] = kep2car([kep ksun]);

%% Velocities
v_in = V_helio_inFb - v_p;

v_out = V_helio_outFb - v_p;

%% Pericenter position and velocities
[~, ~, ~, ~, r_p] = powerGA(V_helio_inFb, V_helio_outFb, FlybyID, ...
    dates(2), 0);

v_p_in = sqrt(norm(v_in)^2+2*mu/r_p);
v_p_out = sqrt(norm(v_out)^2+2*mu/r_p);

%% Hyperbolas
% Saturn SOI data
% Saturn Semimajor axis in AU from https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
a_saturn = 9.53707032;
r_soi_saturn = AU * a_saturn * ( mu / ksun )^( 2/5 );

% Entering hyperbola
e_min = 1 + (r_p * norm(v_in)^2 ) / mu;
delta_min = 2 * asin( 1 / e_min );
DELTA_min = r_p * sqrt( 1 + 2 * ( mu / ( r_p * norm(v_in)^2 ) ) );
theta_inf_min = acos( - 1 / e_min );
beta_min = acos( 1 / e_min );
a_min = DELTA_min / ( e_min^2 - 1 );
b_min = a_min * ( sqrt( e_min^2 - 1 ) );
h_min = sqrt( mu * a_min * ( e_min^2 - 1 ) );

% Exiting hyperbola
e_plus = 1 + ( r_p * norm(v_out)^2 ) / mu;
delta_plus = 2 * asin( 1 / e_plus );
DELTA_plus = r_p * sqrt( 1 + 2 * ( mu / ( r_p * norm(v_out)^2 ) ) );
theta_inf_plus = acos( - 1 / e_plus );
beta_plus = acos( 1 / e_plus );
a_plus = DELTA_plus / ( e_plus^2 - 1 );
h_plus = sqrt( mu * a_plus * ( e_plus^2 - 1 ) );
b_plus = a_plus * ( sqrt( e_plus^2 - 1 ) );

% Velocities at pericenter
vp_min = DELTA_min * norm(v_in) / r_p;
vp_plus = DELTA_plus * norm(v_out) / r_p;

% DeltaV flyby
DELTA_FLYBY = norm(v_out - v_in);
DELTA_VP = abs(vp_plus - vp_min);

% Hyperbola parameters
theta_SOI_min = acos((h_min^2/(mu*r_soi_saturn*e_min))-1/e_min);
theta_SOI_plus = acos((h_plus^2/(mu*r_soi_saturn*e_plus))-1/e_plus);

% Flyby Time
F_min = acosh((cos(theta_SOI_min) + e_min)/(1 + e_min*cos(theta_SOI_min)));
dt_min = sqrt(a_min^3/mu)*(e_min*sinh(F_min)-F_min);
F_plus = acosh((cos(theta_SOI_plus) + e_plus)/(1 + e_plus*cos(theta_SOI_plus)));
dt_plus = sqrt(a_plus^3/mu)*(e_plus*sinh(F_plus)-F_plus);
dt_tot = dt_min+dt_plus;
dt_tot_days = dt_tot/86400;

altitude = r_p - astroConstants(20 + FlybyID);

% Velocity at pericenter vector
v_p_in_vec = vp_min * [0; 1; 0];
v_p_out_vec = vp_plus * [0; 1; 0];

% Integrate the two hyperbola arcs for the 3D plot

[ ~, Y_hyp_in ] = ode113( @(t,y) ode_2bp(t,y,mu), (86400:-1:0), [r_p; 0; 0; 0; vp_min; 0], options );
R_hyp_in = Y_hyp_in(:,1:3);

[ ~, Y_hyp_out ] = ode113( @(t,y) ode_2bp(t,y,mu), (0:1:86400), [r_p; 0; 0; 0; vp_plus; 0], options );
R_hyp_out = Y_hyp_out(:,1:3);

%% Plot

figure;
opt.Position = [0, 0, 0]';
opt.Units = 'km';
opt.FaceAlpha = .75;
opt.RefPlane = 'ecliptic';

ax1 = axes;
h1 = planet3D('Saturn', opt);
hold on
h2 = plot3(ax1, R_hyp_in(:,1), R_hyp_in(:,2), R_hyp_in(:,3),...
    'LineWidth', 2);
hold on
h3 = plot3(ax1, R_hyp_out(:,1), R_hyp_out(:,2), R_hyp_out(:,3),...
    'LineWidth', 2);
h4 = plot3(ax1, R_hyp_out(1,1), R_hyp_out(1,2), R_hyp_out(1,3),...
    'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
hold on
h5 = quiver3(ax1, R_hyp_out(1,1), R_hyp_out(1,2), R_hyp_out(1,3), ...
    (v_p_out_vec(1)-v_p_in_vec(1)), (v_p_out_vec(2)-v_p_in_vec(2)), ...
    (v_p_out_vec(3)-v_p_in_vec(3)), 1e9,'LineWidth', 2, 'color', 'r',...
    'MaxHeadSize', 10);

xlabel(ax1, '$x_{perifocal}$ $[km]$')
ylabel(ax1, '$y_{perifocal}$ $[km]$')
zlabel(ax1, '$z_{perifocal}$ $[km]$')
grid minor
view([0, 90])
legend([h2 h3 h4 h5],'Incoming Hyperbola', 'Outgoing Hyperbola',...
    'Pericenter', '$\Delta V$', 'Location', 'best')

% Display useful data
text = sprintf(['Total time of flight %f hours, pericenter radius %f km,'...
    ' pericenter altitude %f km, SOI radius %f km'], [hours(days(dt_tot_days))...
    , r_p, altitude, r_soi_saturn]);
disp(text)

end
