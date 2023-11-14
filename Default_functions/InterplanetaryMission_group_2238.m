
% This is the code for the Assignment 1. 
% The script gives all the results required:
%   - Porkchop and surface plot of the best time window
%   - Heliocentric transfer trajectory
%   - Planetocentric hyperbolic trajectory

%%
clearvars
close all
clc

%% PATH

addpath('Default_function')
addpath('Default_function\timeConversion\time')
addpath('functions')

%% Plot settings

set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaultlegendInterpreter','latex');
set(0,'defaultAxesFontSize', 11);

%% Initial data

Dep_date_early = datetime([2025 08 02 0 0 0]);
Arr_date_late = datetime([2060 01 29 0 0 0]);
tot_window = Arr_date_late - Dep_date_early;

DepPlanetID = 3;
FbPlanetID = 6;
ArrPlanetID = 61;

%% Loops data

stepDep = 100;
stepFb = 100;
stepArr = 100;

%% Date windows optimization

[dates, minDV, DVArray, datesArray] = DateWindows(Dep_date_early, Arr_date_late,...
    DepPlanetID, FbPlanetID, ArrPlanetID, stepDep, stepFb, stepArr);

%% Plot results 

for windowSelector = 1:size(DVArray, 4)
    for i = 1:size(DVArray, 1)
        for j = 1:size(DVArray, 3)
            DvWindowij(:) = DVArray(i, :, j, windowSelector);
            [minimum, location(i, j, windowSelector)] = min(DvWindowij, [], 'all', 'omitnan');
            if not(isnan(minimum))
                DVArray_plot(i, j, windowSelector) = DVArray(i, location(i, j, windowSelector), j, windowSelector);
                DepArray_plot(i, j, windowSelector) = date2mjd2000(datevec(datesArray(i, location(i, j, windowSelector), j, 1, 1, windowSelector)));
                ArrArray_plot(i, j, windowSelector) = date2mjd2000(datevec(datesArray(i, location(i, j, windowSelector), j, 1, 3, windowSelector)));
            else
                DVArray_plot(i, j, windowSelector) = nan;
                DepArray_plot(i, j, windowSelector) = nan;
                ArrArray_plot(i, j, windowSelector) = nan;
            end
            clear DvWindowij
        end
    end
end

[~, locationWindow] = min(DVArray_plot, [], 'all', 'omitnan');
[~, ~, BestWindow] = ind2sub(size(DVArray_plot), locationWindow);

plotx(:,:) = DepArray_plot(:, :, BestWindow);
ploty(:,:) = ArrArray_plot(:, :, BestWindow);
plotz(:,:) = DVArray_plot(:, :, BestWindow);

[~, pos] = min(plotz, [], 'all', 'omitnan');
[posx, posy] = ind2sub(size(plotz), pos);


surf(plotx, ploty, plotz, 'FaceColor','interp')
hold on
plot3(plotx(posx, posy), ploty(posx, posy), plotz(posx, posy), 'Marker',...
    'o', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
legend('', 'Minimum')
hold off
grid minor
cb2 = colorbar;
datetick('x','yy/mm/dd','keepticks','keeplimits')
datetick('y','yy/mm/dd','keepticks','keeplimits')
set(gca,'XTickLabelRotation',-45)
set(gca,'YTickLabelRotation',45)
xlabel('Departure date [Year/Month/Day]', 'Rotation', 12)
ylabel('Arrival date [Year/Month/Day]','Rotation', -22)
zlabel('Total $\Delta V$ $[km/s]$')

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,-0.003,1]+[0 -.01 0],'Rotation',20)
yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[1,-0.003,1]+[0 -.01 0],'Rotation',-30)

cb2.Label.String = 'Total $\Delta V$ $[km/s]$';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 11;
axis tight
fig = gcf;
set(fig,'PaperOrientation','landscape');
% print(fig, 'ResultFigures/Surf_window3', '-dpdf', '-bestfit')

fig = figure;
ax2 = axes;
ax1 = axes;
ax3 = axes;
contour(ax1, plotx,ploty,(ploty-plotx), 10,'r','ShowText','on', 'LineWidth', 2);
contour(ax2, plotx, ploty, plotz, 50);
plot(ax3, plotx(posx, posy), ploty(posx, posy), 'Marker', 'o',...
    'MarkerFaceColor', 'r', 'MarkerSize', 8)
legend(ax3, 'Minimum')
linkprop([ax1, ax2, ax3],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
ax3.Visible = 'off';
ax3.XTick = [];
ax3.YTick = [];
set([ax1, ax2, ax3], 'Position', [.17 .11 .685 .815]);
cb2 = colorbar(ax2, 'Position', [.88 .11 .0675*0.5 .815]);
datetick(ax2, 'x','yy/mm/dd','keepticks','keeplimits')
datetick(ax2, 'y','yy/mm/dd','keepticks','keeplimits')
set(ax2,'XTickLabelRotation',45)
set(ax2,'YTickLabelRotation',45)
xlabel(ax2, 'Departure date [Year/Month/Day]')
ylabel(ax2, 'Arrival date [Year/Month/Day]')
cb2.Label.String = 'Total $\Delta V$ $[km/s]$';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = 11;
grid minor
fig = gcf;
% print(fig, 'ResultFigures/Porkchop_window3', '-dpdf', '-bestfit')

%% Refinement cycle

[dates_Opt, minDV_Opt] = RefineCycle(dates, minDV, Dep_date_early, Arr_date_late, DepPlanetID, FbPlanetID, ArrPlanetID);

disp(['The minimum delta-v is ', num2str(sum(minDV_Opt)), 'km/s. The optimal departure date is ', string(dates_Opt(1)), ',']);
disp(['the optimal gravity assist date ', string(dates_Opt(2)), ', the optimal arrival date ', string(dates_Opt(3)), '.']);

%% Plot resulting trajectory

plot_helio(dates_Opt, DepPlanetID, FbPlanetID, ArrPlanetID)
fig = gcf;
% print(fig, 'ResultFigures/TransferOrbits', '-dpdf', '-bestfit')

plot_powerGA(dates_Opt, DepPlanetID, FbPlanetID, ArrPlanetID)
fig = gcf;
% print(fig, 'ResultFigures/Hyperbola', '-dpdf', '-bestfit')

%% Save results
% 
% save('Assignment1_results.mat')
% 
