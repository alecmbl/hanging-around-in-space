function [deltaVp, a_in, a_out, delta, r_p] = powerGA(V_helio_in, ...
    V_helio_out, planet_id, GA_date, flagmjd2000)
%
% Evaluation of the main characteristics of a powered fly-by manoeuver
% 
% PROTOTYPE:
%   [deltaVp, a_in, a_out, delta, r_p] = powerGA(V_helio_in, V_helio_out,...
%                                       planet_id, GA_date, flagmjd2000)
%
% DESCRIPTION:
%   The function evaluates the main characteristics for the parabolic orbit
%   around a planet. The fly-by altitude is set at a minimum value of 1e3
%   km to not crash on the planet. This value refers to the Saturn planet's
%   atmosphere. The dates can be either expressed in a vectorial form [Year 
%   Month Days Hours Minutes Seconds], datetime form or modified Julian
%   date. If the latter is the case the flag must be set to 1, otherwise it
%   must be set to 0.
% 
% INPUT:
%   V_helio_in      [3x1]      Heliocentric incoming velocity       [km/s]
%   V_helio_out     [3x1]      Heliocentric outgoing velocity       [km/s]
%   planet_id       [1x1]      Fly-by planet ID                     [-]
%   GA_date         [1x1]      Fly-by date                          [datetime]
%   flagmjd2000     [1x1]      Flag to use Modified Julian days     [-]
%
% OUTPUT:
%   deltaVp         [1x1]   Norm of the powered Delta V required    [km/s]
%   a_in            [1x1]   Incoming hyperbola's semi-major axis    [km]
%   a_out           [1x1]   Outgoing hyperbola's semi-major axis    [km]
%   delta           [1x1]   Deflection angle                        [rad]
%   r_p             [1x1]   Pericenter's radius                     [km]
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

%% Turn off warings of the internal solvers

warning off

%% Date
date_vector = strcmp(class(GA_date), 'double');
date_time = strcmp(class(GA_date), 'datetime');

if flagmjd2000
    
    t_GA = GA_date;

elseif date_vector
    
    [Y, M, D] = ymd(datetime(GA_date));
    [H, m, s] = hms(datetime(GA_date));
    GA_date = [Y M D H m s];

    t_GA = date2mjd2000(GA_date);     % [days]

elseif date_time

    t_GA = date2mjd2000(datevec(GA_date));  % [days]

end


%% Planet

[kep, ~] = uplanet(t_GA, planet_id);
mu = astroConstants(planet_id + 10);
ksun = astroConstants(4);

[~, v_p] = kep2car([kep ksun]);

%% Velocities
v_in = V_helio_in - v_p;
v_out = V_helio_out - v_p;

%%

a_in = mu/norm(v_in)^2;
a_out = mu/norm(v_out)^2;

%% Turning angle

delta = acos(dot(v_in, v_out)/(norm(v_in)*norm(v_out)));

%% R_p

e_minus = @(r_p) 1 + r_p * norm(v_in)^2 / mu;
delta_minus = @(r_p) 2 * asin( 1 / e_minus(r_p) );

e_plus = @(r_p) 1 + r_p * norm(v_out)^2 / mu;
delta_plus = @(r_p) 2 * asin(1/e_plus(r_p));

fun = @(r_p) delta_minus(r_p) / 2 + delta_plus(r_p) / 2 - delta;
options = optimoptions('lsqnonlin','FunctionTolerance',1e-13,'Display','off');

try

    r_p = lsqnonlin(fun, astroConstants(planet_id + 20) + 1e3, 0, inf, options );

catch

    r_p = NaN;

end


if r_p <= (astroConstants(planet_id + 20) + 1e3) || ~isreal(r_p) || isnan(r_p)
    
    deltaVp = nan;
    a_in = nan;
    a_out = nan;
    r_p = nan;
    delta = nan;
    
else

    v_p_in = sqrt(norm(v_in)^2+2*mu/r_p);
    v_p_out = sqrt(norm(v_out)^2+2*mu/r_p);

    deltaVp = abs(v_p_out - v_p_in);

end

end
