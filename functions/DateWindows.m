function [dates, minDV, DVArray, datesArray] = DateWindows(Dep_date_early,...
    Arr_date_late, DepPlanetID, FbPlanetID, ArrPlanetID, stepDep, stepFb,...
    stepArr)
%
% Cycle to find the optimal initial guess of an orbit transfer with a
% powered fly-by in terms of cost of the whole mission
% 
% PROTOTYPE:
%   [dates, minDV, DVArray, datesArray] = DateWindows(Dep_date_early,...
%   Arr_date_late, DepPlanetID, FbPlanetID, ArrPlanetID, stepDep, stepFb,...
%   stepArr)
%
% DESCRIPTION:
%   Finds the best solution in terms of cost of the whole mission for an
%   orbit trasfer with a powered fly-by. It considers the synodic periods
%   to guess the best solution with a coarse grid search. The solution must
%   be refined later to get the actual global minimum for the mission.
% 
% INPUT:
%   Dep_date_early    [1x1]   Earliest Departure date                    [datetime]
%   Arr_date_late     [1x1]   Latest Arrival date                        [datetime]
%   DepPlanetID       [1x1]   Departure planet ID                        [-]
%   FbPlanetID        [1x1]   Fly-by planet ID                           [-]
%   ArrPlanetID       [1x1]   Arrival planet ID                          [-]
%   stepDep           [1x1]   Number of steps for the grid of departure  [-]
%   stepFb            [1x1]   Number of steps for the grid of fly-by     [-]
%   stepArr           [1x1]   Number of steps for the grid of arrival    [-]
% 
% OUTPUT:
%   dates         [3x1]   Best Departure-FlyBy-Arrival dates             [datetime]
%   minDV         [3x1]   Best DeltaV for Departure-Flyby-Arrival        [km/s]
%   DVArray       [121 x 101 x 81 x n] Grid of delta V to construct the
%                                      porkchop plot                     [km/s]
%   datesArray    [121 x 101 x 81 x 1 x 3 x n]  Grid of dates to construct
%                                               the porkchop plot        [datetime]
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



%% Find local minimum first leg inside the first synodic period
[dates1, synodic1, ~, ~, ~] = DateOpt(Dep_date_early, 0, Arr_date_late,...
    stepDep, stepFb, stepArr, DepPlanetID, FbPlanetID, ArrPlanetID,...
    1, 0, 0, 0);

%% Find local minimum second leg inside second synodic period
[dates2, synodic2, ~, ~, ~] = DateOpt(dates1(1,:), 0, Arr_date_late,...
    stepDep, stepFb, stepArr, DepPlanetID, FbPlanetID, ArrPlanetID,...
    2, 0, 0, 0);

%% Find global minimum inside the whole time windows

totdays = Arr_date_late - dates2(1, :);
tot_window = Arr_date_late - Dep_date_early;

Jmax = floor(totdays/synodic2)+1;

dates3 = NaT(3, Jmax);
minDV3 = nan(3, Jmax);
DVArray3 = nan(121, 101, 81, Jmax);
DvTot = nan(Jmax);
datesArray3 = NaT(121, 101, 81, 1, 3, Jmax);

for j = 1:Jmax

    dates(1,:) = dates1(1,:) + (j - 1) * synodic1 * ceil(synodic2/synodic1);
    dates(2,:) = dates2(2,:) + (j - 1) * synodic2;

    if dates(1,:) <= Arr_date_late && dates(2,:) <= Arr_date_late

    [dates3(:, j), ~, minDV3(:, j), DVArray3(:, :, :, j),...
        datesArray3(:, :, :, :, :, j)] = DateOpt(dates(1,:), dates(2,:),...
        Arr_date_late, 120, 100, 80, DepPlanetID, FbPlanetID,...
        ArrPlanetID, 3, synodic2/100, synodic1/100, tot_window/20);
    DvTot(j) = sum(minDV3(:, j));

    end

end

%% Find the absolute minimum

[~,location] = min(DvTot, [],'all', 'omitnan');
dates4 = dates3(:, location);
minDV4 = minDV3(:, location);

%% Output results

dates = dates4;
minDV = minDV4;
DVArray = DVArray3;
datesArray = datesArray3;

end