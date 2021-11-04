function [AQMData, risingData, risingModel, fallingData, fallingModel] = AQM_plotter()
%% Function to import and plot air quality monitor (AQM) data files.
% All files are expected to be 5 columns of floating points, containing the
% time, temperature, relative humidity, VOCs concentration and particle
% count. The function assumes the sampling interval is 1 seconds, i.e. data
% were saved to the PC rather than the memory card.

%% My generic preferred plotting options.
addpath(cd);
set(groot, 'defaultFigureWindowStyle', 'Docked','defaultFigureColor',...
    'White', 'defaultAxesFontSize', 20, 'defaultAxesFontName', 'Arial',...
    'defaultLineLineWidth', 1.3, 'defaultAxesTickDir', 'Out',...
    'defaultAxesTickDirMode', 'Manual', 'defaultAxesLineWidth', 2,...
    'defaultFigureColormap', gray(64));

%% Select and import files.
[files, path] = uigetfile('*.txt', 'MultiSelect', 'On'); % Choose files.
cd(path); % Set current folder to data folder.
if ischar(files) == 1 % If a single file is selected.
    files = {files}; % Convert it to a single cell
end
AQMData = zeros(1, 5); % Matrix to receive data.
warning off
for iFiles = 1 : numel(files) % Iterate through files.
    disp(files{iFiles}); % State file loaded.
    fileData = table2array(readtable(strcat(path, files{iFiles}))); % Read in data.
    fileData(:, 1) = fileData(:, 1) + AQMData(end, 1); % Make time consecutive.
    AQMData = [AQMData; fileData]; % Append to full data matrix.
end
AQMData(1, :) = []; % Remove starting values.
AQMData(:, 1) = (AQMData(:, 1)) ./ 60; % Convert time to mintues.
AQMData(isnan(AQMData(:, 4)), 4) = 50; % Set saturated VOC data to 50;

%% Plot data.
figure;
hold on
box off
plot(AQMData(:, 1), AQMData(:, 4), 'Color', 'k'); % Plot VOC concentration.
ylabel('VOCs [ppm]');
yyaxis right
plot(AQMData(:, 1), AQMData(:, 5)); % Plot particle count.
ylabel('Particle count [\mug/m^{3}]');
xlabel('Time [min]');

%% Find indices of VOC peak feature.
ratData = 1 ./ AQMData(:, 4); % Get ratio data to find peaks with.
ratData(ratData == Inf) = 0; % Remove infinite values.
figure;
subplot(2, 2, 1);
findpeaks(ratData, 'NPeaks', 1, 'MinPeakProminence', 0.01); % Find onset.
[~, VRiseMin] = findpeaks(ratData, 'NPeaks', 1, 'MinPeakProminence', 0.01);
VRiseMin = VRiseMin - 1;
subplot(2, 2, 2);
findpeaks(AQMData(:, 4), 'NPeaks', 1); % Find rising peak.
[~, VRiseMax] = findpeaks(AQMData(:, 4), 'NPeaks', 1);
subplot(2, 2, 3);
findpeaks(flipud(AQMData(:, 4)), 'NPeaks', 1,...
    'MinPeakProminence', 50); % Find falling peak.
[~, VFallMax] = findpeaks(flipud(AQMData(:, 4)), 'NPeaks', 1,...
    'MinPeakProminence', 50);
VFallMax = size(AQMData, 1) - VFallMax; % Flip index to correct.
subplot(2, 2, 4);
findpeaks(flipud(ratData), 'NPeaks', 1, 'MinPeakProminence', 0.5); % Find end.
[~, VFallMin] = findpeaks(flipud(ratData), 'NPeaks', 1, 'MinPeakProminence', 0.5);
VFallMin = size(AQMData, 1) - (VFallMin - 1); % Flip index to correct.
risingData = [AQMData(VRiseMin : VRiseMax, 1), AQMData(VRiseMin : VRiseMax, 4)];
fallingData = [AQMData(VFallMax : VFallMin, 1), AQMData(VFallMax : VFallMin, 4)];

%% Fit exponentials to VOC trace to estimate peak value.
fittingExt = risingData(1, 1) : 0.005 : fallingData(end, 1); % To extrapolate data.
risingModel = fit(risingData(:, 1), risingData(:, 2), 'exp1'); % Rising model.
risingFit = risingModel(fittingExt); % Fit rising data.
fallingModel = fit(fallingData(:, 1), fallingData(:, 2), 'exp2'); % Falling model.
fallingFit = fallingModel(fittingExt); % Fit falling data.
[~, extCross] = min(abs(risingFit - fallingFit)); % Find intersection.
extCross = risingFit(extCross); % Find crossing value.
figure;
plot(risingData(:, 1), risingData(:, 2), 'o', 'Color', 'k');
hold on
box off
plot(fittingExt, risingFit, 'Color', 'k');
plot(fallingData(:, 1), fallingData(:, 2), 'o', 'Color', 'r');
plot(fittingExt, fallingFit, 'Color', 'b');
ylim([0, extCross + 20]);
ylabel('VOCs [ppm]');
xlabel('Time [min]');
legend('Rising data', 'Rising fit', 'Falling data', 'Falling fit');
legend box off
title(strcat('Intersection at', {' '}, num2str(ceil(extCross)), ' ppm.'));

end