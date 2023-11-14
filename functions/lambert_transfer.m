function [DeltaV1, DeltaV2, DeltaV_tot, V_helio_in, V_helio_out] = ...
    lambert_transfer(t1_d, t2_d, ToF, planet1_id, planet2_id, mu_p,...
    flagmjd2000)
%
%   Cost of a transfer between two planets with a lambert's arc
% 
% PROTOTYPE:
%   [DeltaV1, DeltaV2, DeltaV_tot, V_helio_in, V_helio_out] = 
%        lambert_transfer(t1_d, t2_d, ToF, planet1_id, planet2_id, mu_p,...
%        flagmjd2000)
%
% DESCRIPTION:
%   The function evaluates the cost for an orbit transfer between two
%   planets in the solar system using a lambert solver. The dates can be 
%   either expressed in a vectorial form [Year Month Days Hours Minutes 
%   Seconds], datetime form or modified Julian date. If the latter is the 
%   case the flag must be set to 1, otherwise it must be set to 0.
% 
% INPUT:
%   t1_d            [1x1]      Departure date                           [datetime]
%   t2_d            [1x1]      Arrival date                             [datetime]
%   ToF             [1x1]      Time of flight                           [s]
%   planet1_id      [1x1]      Departure planet ID                      [-]
%   planet2_id      [1x1]      Arrival planet ID                        [-]
%   mu_p            [1x1]      Gravitational constant of the primary    [km^2/s^3]
%   flagmjd2000     [1x1]      Flag to use Modified Julian days         [-]
%
% OUTPUT:
%   DeltaV1         [1x1]   Norm of the Delta V required at departure   [km/s]
%   DeltaV2         [1x1]   Norm of the Delta V required at arrival     [km/s]
%   V_helio_in      [3x1]   Heliocentric velocity required at departure [km/s]
%   V_helio_out     [3x1]   Heliocentric velocity required at arrival   [km/s]
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



% This function computes the lambert solution for two given points and a
% time of flight and plots the resulting transfer arc

%% Date conversions
date_vector = strcmp(class(t1_d), 'double');
date_time = strcmp(class(t1_d), 'datetime');

if date_time

    t1mjd2000 = date2mjd2000(datevec(t1_d));    % [days]
    t2mjd2000 = date2mjd2000(datevec(t2_d));    % [days]

elseif flagmjd2000
    
    t1mjd2000 = t1_d;
    t2mjd2000 = t2_d;

elseif date_vector
    
    [Y1, M1, D1] = ymd(datetime(t1_d));
    [H1, m1, s1] = hms(datetime(t1_d));
    t1_d = [Y1 M1 D1 H1 m1 s1];

    [Y2, M2, D2] = ymd(datetime(t2_d));
    [H2, m2, s2] = hms(datetime(t2_d));
    t2_d = [Y2 M2 D2 H2 m2 s2];

    t1mjd2000 = date2mjd2000(t1_d);             % [days]
    t2mjd2000 = date2mjd2000(t2_d);             % [days]

end

%% Boundary conditions
% Ephemerides

[kep1, ~] = uplanet(t1mjd2000, planet1_id);

if planet2_id < 11
    [kep2, ~] = uplanet(t2mjd2000, planet2_id);
else
    [kep2, ~, ~] = ephNEO(t2mjd2000, planet2_id);
end


% Initial and final position 
[r1, v1] = kep2car([kep1 mu_p]);
[r2, v2] = kep2car([kep2 mu_p]);

%% Transfer orbit

orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 1;

[~,~,~,~,v1_transf, v2_transf,~,~] = lambertMR(r1, r2, ToF, mu_p, ...
    orbitType, Nrev, Ncase, optionsLMR);


%% DeltaV

DeltaV(1) = abs(norm(v1_transf' - v1));
DeltaV1 = DeltaV(1);

DeltaV(2)= abs(norm(v2_transf' - v2));
DeltaV2 = DeltaV(2);

DeltaV_tot = sum(DeltaV);

%% Velocities

V_helio_in = v1_transf';
V_helio_out = v2_transf';

end