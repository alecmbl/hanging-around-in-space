function [dates, synodic, minDV, DVarray, datesArray] = ...
    DateOpt(Dep_date_early, Fb_date, Arr_date_late, stepDep, nstepToF1,...
    nstepToF2, DepartureID, FlybyID, ArrivalID, flagSyn, stepFb,...
    stepDep_time, stepArr)

%
% Triple cycle to find optimal mission cost and corresponding dates for an
% orbit transfer with powered fly-by
%
% PROTOTYPE:
%   [dates, synodic, minDV, DVarray, datesArray] = DateOpt(Dep_date_early,...
%                              Fb_date, Arr_date_late, stepDep, nstepToF1,...
%                              nstepToF2, DepartureID, FlybyID, ArrivalID,...
%                              flagSyn, stepFb, stepDep_time, stepArr)
%
% DESCRIPTION:
%   Finds the best solution for an orbit transfer with powered fly-by using
%   a grid search. The flag must be set to 1, 2, 3 or 4. 1 = best solution
%   for the first leg, 2 = best solution for the second leg, 3 = best
%   overall solution, 4 = best solution for the first leg with a time
%   window set by the inputs. The dates can be either expressed in a
%   vectorial form [Year Month Days Hours Minutes Seconds] or datetime form.
%   The earliest and latest dates are set by default in the triple nested
%   loop as [2025 08 02 0 0 0] and [2060 01 29 0 0 0]. Any non-feasible
%   solution is replaced with NaN or NaT in the arrays and matrices.
%
% INPUT:
%   Dep_date_early    [1x1]   Earliest Departure date                    [datetime]
%   Fb_date           [1x1]   Fly-by date                                [datetime]
%   Arr_date_late     [1x1]   Latest Arrival date                        [datetime]
%   stepDep           [1x1]   Number of steps for the grid of departure  [-]
%   nstepToF1         [1x1]   Number of steps for the grid of fly-by     [-]
%   nstepToF2         [1x1]   Number of steps for the grid of arrival    [-]
%   DepartureID       [1x1]   Departure planet ID                        [-]
%   FlybyID           [1x1]   Fly-by planet ID                           [-]
%   ArrivalID         [1x1]   Arrival planet ID                          [-]
%   flagSyn           [-]     Flag to choose the cases listed in
%                              description                               [-]
%   stepFb            [1x1]   Time lenght of the step in the grid for
%                               fly-by dates                             [duration]
%   stepDep_time      [1x1]   Time lenght of the step in the grid for
%                               departure dates                          [duration]
%   stepArr           [1x1]   Time lenght of the step in the grid for
%                               arrival dates                            [duration]
%
% OUTPUT:
%   dates         [3x1]   Best Departure-FlyBy-Arrival dates             [datetime]
%   synodic       [1x1]   Synodic period used in computation if any      [duration]
%   minDV         [3x1]   Best DeltaV for Departure-Flyby-Arrival        [km/s]
%   DVArray       [(stepDep+1) x (nstepToF1+1) x (nstepToF1+1) x n] Grid
%                                       of delta V                       [km/s]
%   datesArray    [(stepDep+1) x (nstepToF1+1) x (nstepToF1+1) x n] Grid
%                                       of dates                         [datetime]
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


%% Dates manipulation
date_vector = strcmp(class(Dep_date_early), 'double');

if date_vector
    Dep_date_early = datetime(Dep_date_early);
    Arr_date_late = datetime(Arr_date_late);
    Fb_date = datetime(Fb_date);
end

Dep_date_early_vec = datevec(Dep_date_early);       % vectorial
Arr_date_late_vec = datevec(Arr_date_late);         % vectorial

%% Planet sidereal and synodic periods
% Initial positions
DepMJD2000 = date2mjd2000(Dep_date_early_vec);
[kepDep, ksun] = uplanet(DepMJD2000, DepartureID);
[kepFlyBy, ~] = uplanet(DepMJD2000, FlybyID);
[kepArr, ~, ~] = ephNEO(DepMJD2000, ArrivalID);

% Sidereal periods
DeparturePeriod = 2 * pi * sqrt(kepDep(1)^3 / ksun);
FlybyPeriod = 2 * pi * sqrt(kepFlyBy(1)^3 / ksun);
ArrivalPeriod = 2 * pi * sqrt(kepArr(1)^3 / ksun);

% Synodic periods
SynPeriodDep_FB = seconds((FlybyPeriod * DeparturePeriod)/(abs(FlybyPeriod - DeparturePeriod)));
SynPeriodFB_Arr = seconds((ArrivalPeriod * FlybyPeriod)/(abs(ArrivalPeriod - FlybyPeriod)));

%% Departure date interval

if flagSyn == 0

    [synodic, index] = min(SynPeriodDep_FB, SynPeriodFB_Arr);
    flagSyn = index;

end

if flagSyn == 1

    synodic = SynPeriodDep_FB;
    DepArray = linspace(Dep_date_early, Dep_date_early + synodic, stepDep + 1);   % [datetime]
    FbArray = linspace(Dep_date_early, Arr_date_late, nstepToF1 + 1);             % [datetime]
    dimDep = length(DepArray);
    dimFb = length(FbArray);
    ArrArray = linspace(Dep_date_early, Arr_date_late, nstepToF2 + 1);             % [datetime]

elseif flagSyn == 2

    synodic = SynPeriodFB_Arr;
    DepArray = Dep_date_early;
    FbArray = linspace(Dep_date_early, Dep_date_early + synodic, nstepToF1 + 1);    % [datetime]
    dimDep = 1;
    dimFb = length(FbArray);
    ArrArray = linspace(Dep_date_early, Arr_date_late, nstepToF2 + 1);             % [datetime]

elseif flagSyn == 3 || flagSyn == 4

    synodic = SynPeriodFB_Arr;
    DepArray = linspace(Dep_date_early - stepDep/2*stepDep_time, Dep_date_early + stepDep/2*stepDep_time, stepDep + 1);    % [datetime]
    FbArray = linspace(Fb_date - nstepToF1/2*stepFb, Fb_date + nstepToF1/2*stepFb, nstepToF1 + 1);    % [datetime]
    ArrArray = linspace(Arr_date_late - nstepToF2/2*stepArr, Arr_date_late + nstepToF2/2*stepArr, nstepToF2 + 1);    % [datetime]
    dimDep = length(DepArray);
    dimFb = length(FbArray);

end

dimArr = length(ArrArray);

%% Preallocating arrays
ToF1 = nan(dimDep, dimFb, dimArr);
ToF2 = nan(dimDep, dimFb, dimArr);
DeltaVLeg1 = nan(dimDep, dimFb, dimArr);
DeltaVLeg2 = nan(dimDep, dimFb, dimArr);
DeltaVFB = nan(dimDep, dimFb, dimArr);
DeltaVman = nan(dimDep, dimFb, dimArr, 3);
DeltaV_tot = nan(dimDep, dimFb, dimArr);
time = NaT(dimDep, dimFb, dimArr, 1, 3);

%% Triple nested loop

parfor i = 1:dimDep

    Dep_date = DepArray(i);                 % [datetime]

    if Dep_date >= datetime([2025 08 02 0 0 0])

        for j = 1:dimFb

            Fb_date = FbArray(j);               % [datetime]

            if Fb_date > Dep_date && Fb_date < datetime([2060 01 29 0 0 0])

                for k = 1:dimArr

                    Arr_date = ArrArray(k);             % [datetime]

                    if Arr_date > Fb_date && Arr_date <= datetime([2060 01 29 0 0 0])

                        time(i, j, k, :, :) = [Dep_date Fb_date Arr_date];

                        ToF1(i, j, k) = seconds(Fb_date - Dep_date);
                        ToF2(i, j, k) = seconds(Arr_date - Fb_date);

                        if ToF1(i, j, k) > 0 && ToF2(i, j, k) > 0

                            [DeltaVLeg1(i, j, k), ~, ~, ~, V_helio_in] = lambert_transfer(Dep_date, Fb_date, ToF1(i, j, k), DepartureID, FlybyID, ksun, 0);
                            [~, DeltaVLeg2(i, j, k), ~, V_helio_out, ~] = lambert_transfer(Fb_date, Arr_date, ToF2(i, j, k), FlybyID, ArrivalID, ksun, 0);
                            [DeltaVFB(i, j, k), ~, ~, ~, ~] = powerGA(V_helio_in, V_helio_out, FlybyID, Fb_date, 0);

                            if ~isnan(DeltaVLeg1(i, j, k)) && ~isnan(DeltaVLeg2(i, j, k)) && ~isnan(DeltaVFB(i, j, k))

                                DeltaVman(i, j, k, :) = [DeltaVLeg1(i, j, k)  DeltaVFB(i, j, k) DeltaVLeg2(i, j, k)  ];
                                DeltaV_tot(i, j, k) = sum(DeltaVman(i, j, k, :));

                            end

                        else

                            DeltaV_tot(i, j, k) = NaN;

                        end

                    end

                end

            end

        end

    end
end

%% Best solution research

if flagSyn == 1 || flagSyn == 4

    % Leg 1 minimum

    DeltaVLeg1(:,:,:) = DeltaVman(:,:,:,1);
    [~, locMinLeg1] = min(DeltaVLeg1, [], 'all', 'omitnan');
    [idx1minLeg1,idx2minLeg1, idx3minLeg1] = ind2sub(size(DeltaVLeg1), locMinLeg1);
    DepOpt1 = time(idx1minLeg1, idx2minLeg1, idx3minLeg1, :, 1);
    FbOpt1 = time(idx1minLeg1, idx2minLeg1, idx3minLeg1, :, 2);
    ArrOpt1 = time(idx1minLeg1, idx2minLeg1, idx3minLeg1, :, 3);
    dates = [DepOpt1; FbOpt1; ArrOpt1];
    minDV(:, :, :) = DeltaVman(idx1minLeg1, idx2minLeg1, idx3minLeg1, :);

elseif flagSyn == 2

    % Leg 2 minimum

    DeltaVLeg2(:,:,:) = DeltaVman(:,:,:,3);
    [~, locMinLeg2] = min(DeltaVLeg2, [], 'all', 'omitnan');
    [idx1minLeg2,idx2minLeg2, idx3minLeg2] = ind2sub(size(DeltaVLeg2), locMinLeg2);
    DepOpt2 = time(idx1minLeg2, idx2minLeg2, idx3minLeg2, :, 1);
    FbOpt2 = time(idx1minLeg2, idx2minLeg2, idx3minLeg2, :, 2);
    ArrOpt2 = time(idx1minLeg2, idx2minLeg2, idx3minLeg2, :, 3);
    dates = [DepOpt2; FbOpt2; ArrOpt2];
    minDV(:, :, :) = DeltaVman(idx1minLeg2, idx2minLeg2, idx3minLeg2, :);

elseif flagSyn == 3

    % Overall minimum

    [~,location] = min(DeltaV_tot(:), [], 'omitnan');
    [idx1,idx2,idx3] = ind2sub(size(DeltaV_tot),location);
    DepOpt3 = time(idx1, idx2, idx3, :, 1);
    FbOpt3 = time(idx1, idx2, idx3, :, 2);
    ArrOpt3 = time(idx1, idx2, idx3, :, 3);
    dates = [DepOpt3; FbOpt3; ArrOpt3];
    minDV(:, :, :) = DeltaVman(idx1, idx2, idx3, :);

end

DVarray = DeltaV_tot;
datesArray = time;
synodic = years(years(synodic));


end

