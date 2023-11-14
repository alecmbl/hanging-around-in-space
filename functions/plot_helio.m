function plot_helio(dates, DepartureID, FlybyID, ArrivalID)
%
% Plot of the transfer trajectory of an orbit transfer with powered fly-by
% 
% PROTOTYPE:
%   plot_helio(dates, DepartureID, FlybyID, ArrivalID)
%
% DESCRIPTION:
%   Plots the transfer trajectory, the planet's position in the dates
%   specified in the input and the planet's orbits with a dashed line.
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
elseif date_time
    Dep_date = dates(1);
    Arr_date = dates(3);
    Fb_date = dates(2);
end

%% Constants

AU = astroConstants(2);

%% Planets orbits

% Initial positions
DepMJD2000 = date2mjd2000(datevec(Dep_date));
FlyByMJD2000 = date2mjd2000(datevec(Fb_date));
ArrMJD2000 = date2mjd2000(datevec(Arr_date));

[kepDep, ksun] = uplanet(DepMJD2000, DepartureID);
[kepFlyBy, ~] = uplanet(FlyByMJD2000, FlybyID);
[kepArr, ~, ~] = ephNEO(ArrMJD2000, ArrivalID);

[r_Dep, v_Dep] = kep2car([kepDep ksun]);
[r_FlyBy, v_FlyBy] = kep2car([kepFlyBy ksun]);
[r_Arr, v_Arr] = kep2car([kepArr ksun]);

% Propagate orbits
[RR_Dep, ~] = two_body_problem_fun( r_Dep, v_Dep, ksun);
RR_Dep = RR_Dep./AU;

[RR_FlyBy, ~] = two_body_problem_fun( r_FlyBy, v_FlyBy, ksun);
RR_FlyBy = RR_FlyBy./AU;

[RR_Arr, ~] = two_body_problem_fun( r_Arr, v_Arr, ksun);
RR_Arr = RR_Arr./AU;

%% Transfer trajectories

% Transfer arc 1
[~, ~, ~, V_helio_inDep, ~] = lambert_transfer(Dep_date, Fb_date,...
    seconds(Fb_date - Dep_date), DepartureID, FlybyID, ksun, 0);

tspan = linspace( 0, seconds(Fb_date - Dep_date), 1000 ); % Time span

% Perform the integration
[ ~, Y_Leg1 ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, [r_Dep, V_helio_inDep], options );
R_Leg1 = Y_Leg1(:,1:3)./AU;


% Transfer arc 2
[~, ~, ~, V_helio_outFb, ~] = lambert_transfer(Fb_date, Arr_date,...
    seconds(Arr_date - Fb_date), FlybyID, ArrivalID, ksun, 0);

tspan = linspace( 0, seconds(Arr_date - Fb_date), 1000 ); % Time span
[ ~, Y_Leg2 ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, [r_FlyBy, V_helio_outFb], options );
R_Leg2 = Y_Leg2(:,1:3)./AU;



% Plot

figure;
opt1.Position = [R_Leg1(1,1); R_Leg1(1,2); R_Leg1(1,3)];
opt1.Units = 'AU';
opt1.Clipping = 'on';

opt2 = opt1;
opt2.Position = [R_Leg2(1,1), R_Leg2(1,2), R_Leg2(1,3)]';

opt3 = opt1;
opt3.Position = [0; 0; 0];

ax = axes;


h1 = plot3(ax, RR_Dep(:,1), RR_Dep(:,2), RR_Dep(:,3), '--', 'LineWidth', 1);
hold on
h2 = plot3(ax, RR_FlyBy(:,1), RR_FlyBy(:,2), RR_FlyBy(:,3), '--', ...
    'LineWidth', 1);
h3 = plot3(ax, RR_Arr(:,1), RR_Arr(:,2), RR_Arr(:,3), '--', 'LineWidth', 1);

h4 = plot3(ax, R_Leg1(:,1), R_Leg1(:,2), R_Leg1(:,3), 'LineWidth', 3);
h5 = plot3(ax, R_Leg2(:,1), R_Leg2(:,2), R_Leg2(:,3), 'LineWidth', 3);


hold on
h6 = plot3(ax, R_Leg1(1,1), R_Leg1(1,2), R_Leg1(1,3), 'o', ...
    'MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
h7 = plot3(R_Leg2(1,1), R_Leg2(1,2), R_Leg2(1,3), 'o', 'MarkerSize',...
    5, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
h8 = plot3(ax, R_Leg2(end,1), R_Leg2(end,2), R_Leg2(end,3), 'o',...
    'MarkerSize', 5, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');

Earth = planet3D('Earth', opt1);
Saturn = planet3D('Saturn', opt2);
Sun = planet3D('Sun', opt3);

grid minor
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
legend([h1, h2, h3, h4, h5, h6, h7, h8], 'Earth orbit', 'Saturn orbit', 'NEO 61 orbit',...
    'First leg transfer', 'Second leg transfer', 'Departure', 'Fly-by',...
    'Arrival', 'Location', 'northeast')
view([0, 90])
axis auto

end