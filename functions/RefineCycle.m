function [dates, minDV] = RefineCycle(dates_old, minDV_old, ...
    Dep_date_early, Arr_date_late, DepartureID, FlybyID, ArrivalID)
%
% Refinement cycle to get minimum Delta V and corresponding dates
% 
% PROTOTYPE:
%   [dates, minDV] = RefineCycle(dates_old, minDV_old, Dep_date_early,...
%                   Arr_date_late, DepartureID, FlybyID, ArrivalID)
%
% DESCRIPTION:
%   Finds the best solution with a relative tollerance of 1e-5 for an
%   orbit trasfer with a powered flyby in a restricted time window. The
%   amount of steps considered for the grid search is higher for inner
%   planets and lower for outer planets in the solar system frame of
%   reference. The chosen values refer to the system Earth-Saturn-NEO61 in
%   which the flyby planet is the outermost and the Earth is the innermost.
% 
% INPUT:
%   dates_old         [3x1]   Initial Departure-FlyBy-Arrival dates      [datetime]
%   minDV_old         [3x1]   Initial DeltaV for Departure-Flyby-Arrival [km/s]
%   Dep_date_early    [1x1]   Earliest Departure date                    [datetime]
%   Arr_date_late     [1x1]   Latest Arrival date                        [datetime]
%   DepartureID       [1x1]   Departure planet ID                        [-]
%   FlybyID           [1x1]   Fly-by planet ID                           [-]
%   ArrivalID         [1x1]   Arrival planet ID                          [-]
% 
% OUTPUT:
%   dates         [3x1]   Best Departure-FlyBy-Arrival dates             [datetime]
%   minDV         [3x1]   Best DeltaV for Departure-Flyby-Arrival        [km/s]
%
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

%% Time window

tot_window = days(days(Arr_date_late - Dep_date_early));

%% Planet sidereal and synodic periods
% Initial positions
DepMJD2000 = date2mjd2000(datevec(Dep_date_early));
[kepDep, ksun] = uplanet(DepMJD2000, DepartureID);
[kepFlyBy, ~] = uplanet(DepMJD2000, FlybyID);
[kepArr, ~, ~] = ephNEO(DepMJD2000, ArrivalID);

% Sidereal periods
DeparturePeriod = 2 * pi * sqrt(kepDep(1)^3 / ksun);
FlybyPeriod = 2 * pi * sqrt(kepFlyBy(1)^3 / ksun);
ArrivalPeriod = 2 * pi * sqrt(kepArr(1)^3 / ksun);

% Synodic periods
synodic1 = seconds((FlybyPeriod * DeparturePeriod)/(abs(FlybyPeriod - DeparturePeriod)));
synodic2 = seconds((ArrivalPeriod * FlybyPeriod)/(abs(ArrivalPeriod - FlybyPeriod)));

counter = 1;
stepDep = [];
nstepToF1 = [];
stepDep(1) = 40;
nstepToF1(1) = 40;
nstepToF2(1) = 40*4;
error = [];
error(1) = 1e3;
counter_error = 0;
nstep = 20;

while error(end) > 1e-5
    [dates_new(:,counter), ~, minDV_new(:,counter)] = DateOpt(dates_old(1,:), dates_old(2, :), dates_old(3, :),...
        4*nstep, nstep, 2*nstep, 3, 6, 61, 3, synodic2/nstepToF1(end), synodic1/stepDep(end), tot_window/(nstepToF2(end)*14));

    error(counter) = abs(sum(minDV_new(:,end)) - sum(minDV_old));

    if error(end) == 0


        if counter_error <= 0

            counter_error = 1;

        else

            counter_error = counter_error + 1;

        end

        error(end) = 1e3;

        if counter > 1

            stepDep(end + 1) = stepDep(end) * 2;
            nstepToF1(end + 1) = nstepToF1(end) * 2;
            nstepToF2(end + 1) = nstepToF2(end) * 2;

        end

    else

        counter_error = 0;

    end

    dates_old = dates_new(:,end);
    minDV_old = minDV_new(:,end);

    counter = counter + 1;

    if counter_error == 30
        break
    end

    if counter == 100
        break
    end

end

dates = dates_old;
minDV = minDV_old;

end
