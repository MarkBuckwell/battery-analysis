
%% Mark Buckwell, Electrochemical Innovation Lab, University College London
% THT_thermal_analysis, v1.0 2021

%% Sort out the starting workspace.
clear
clc
close all
fclose('all');
addpath(cd);

%% My generic preferred plotting options. 
set(groot, 'defaultFigureWindowStyle', 'Docked','defaultFigureColor', 'White',...
    'defaultAxesFontSize', 16, 'defaultAxesFontName', 'Arial',...
    'defaultLineLineWidth', 1.3, 'defaultAxesTickDir', 'Out',...
    'defaultAxesTickDirMode', 'Manual', 'defaultAxesLineWidth', 2,...
    'defaultFigureColormap', gray(64));

%% This script will process, analyse and plot THT ARC data.

%% Initialise parameters and counters.
headers = {'time', 'cap', 'crate', 'pressure', 'prate', 'top', 'side',...
    'bottom', 'body', 'base', 'rPower', 'tPower', 'sPower', 'bPower', 'mode'};
startTemp = 45; % Offset time to zero at this temperature.
ventStart = 75; % Ignore this many samples at the start to ignore misc peaks.
ventPeak = 2; % Height of vent peak in temp rate.
ventOffset = 3; % Number of samples to offset back from vent peak.
TRPeak = 100; % Height of thermal runaway peak in temp rate.
TROffset = 1; % Number of samples to offset back from TR peak.
modelType = 'poly2';
disThresh = 1500; % Threshold TC value to consider disconnectied.

%% Choose folder of files to work through.
[fName, fLoc] = uigetfile('*.DAT', 'Multiselect', 'Off'); % Choose file.
cd(fLoc);

%% Open file and process data.
fChoice = fopen(fName); % Open file.
fgetl(fChoice); % Ignore first line.
nCols = textscan(fgetl(fChoice), '%f', 'Delimiter', ' ');
nCols = numel(nCols{1}); % Get number of columns.
frewind(fChoice);
data = array2table(cell2mat(textscan(fChoice, repmat('%f', 1, nCols),...
    'Delimiter', ' ', 'HeaderLines', 1)), 'VariableNames', headers); % Get data.
data.crate(data.cap > disThresh) = NaN; % Remove cap rate TC disconnection artefacts.
data.cap(data.cap > disThresh) = NaN; % Remove cap TC disconnection artefacts.
data.body(data.body > 1500); % Remove body TC disconnection datapoints.
data.base(data.base > 1500); % Remove base TC disconnection datapoints.
[~, tPeak] = max(data.cap); % Find peak temperature index.
[~, tOffset] = min(abs(data.cap(1 : tPeak) - startTemp)); % Find start temp.
data.time = data.time - data.time(tOffset); % Offset time.
[~, TRInd] = findpeaks(data.crate(1 : tPeak),...
    'NPeaks', 1, 'MinPeakHeight', TRPeak); % Get thermal runaway index.
TRInd = TRInd - TROffset; % Offset thermal runaway index.
TRTime = data.time(TRInd); % Get thermal runaway time.
TRTemp = data.cap(TRInd); % Get thermal runaway temperature.
TRError = abs(data.cap(TRInd)...
    - data.cap(TRInd -1)) / 2; % Estimate TR temp. error.
ventRate = data.crate; % Get vent profile from cap data.
figure;
subplot(2, 1, 2);
findpeaks(data.crate(1 : tPeak), data.time(1 : tPeak),...
    'NPeaks', 1, 'MinPeakHeight', TRPeak);
box off
xlabel('Time [min]');
ylabel('d^{2}T/dt^{2}');
iVent = 0; % Initialise venting check.
disp('Checking cap data.');
ventLegend = {'Cap'};
subplot(2, 1, 1);
box off
xlabel('Time [min]');
ylabel('dT/dt');
hold on
while iVent < 3
    try % First go looking for the venting in the cap profile.
        
        %% Plot peak finding.
        findpeaks(-1 .* gradient(ventRate(ventStart : tPeak)),...
            data.time(ventStart : tPeak), 'NPeaks', 1, 'MinPeakHeight', ventPeak);
        [~, ventInd] = findpeaks(-1 .* gradient(ventRate(ventStart : tPeak)),...
            'NPeaks', 1, 'MinPeakHeight', ventPeak); % Get vent index
        
        %% Get event times.
        ventInd = ventInd - ventOffset + ventStart; % Offset venting index.
        ventTime = data.time(ventInd); % Get vent time.
        ventTemp = data.cap(ventInd); % Get vent temperature.
        ventError = abs(data.cap(ventInd + 1)...
            - data.cap(ventInd - 1)) / 2; % Estimate vent temp. error.
                
        %% Fit an exponential to the region between venting and thermal runaway.
        tExp = data.time(ventInd + 1 : TRInd - 1);
        TExp = data.cap(ventInd + 1 : TRInd - 1);
        rExp = data.crate(ventInd + 1 : TRInd - 1);
        accModel = fit(tExp, TExp, modelType);
        rateModel = fit(tExp, rExp, modelType);

        disp('Finished checking for venting.');
        iVent = 3; % Stop looking if no errors.
    catch % If no venting is found in the cap data.
        disp('No vent detected.');
        if iVent == 0 % Use body data instead.
            ventRate = gradient(data.body) ./ gradient(data.time);
            iVent = 1; % Increment vent check.
            disp('Checking body data.');
            ventLegend{2} = 'Body';
        elseif iVent == 1 % Finally try base data.
            ventRate = gradient(data.base) ./ gradient(data.time);
            iVent = 2; % Will end vent checking on next loop.
            disp('Checking base data.');
            ventLegend{3} = 'Base';
        elseif iVent == 2
            iVent = 3; % End vent checking.
            disp('Finished checking for venting.');
        end
    end
end
ylim([-10, 10]);
legend(ventLegend, 'Location', 'NorthWest');
legend box off

%% Plot section between venting and thermal runaway.
figure;
hold on
yyaxis left
plot(tExp, TExp, 'o', 'Color', 'bl');
plot(tExp, accModel(tExp), '-', 'Color', 'bl');
set(gca, 'YScale', 'log');
xlabel('Time [min]');
ylabel('Temperature [˚C]');
yyaxis right
plot(tExp, rExp, 'o', 'Color', 'r');
plot(tExp, gradient(rateModel(tExp))./gradient(tExp), '-', 'Color', 'r');
ylabel('Heating rate [˚C/min]');
set(gca, 'YScale', 'log');

%% Plot temperature profiles and located events.
figure;
subplot(2, 1, 1);
hold on
plot(data.time, data.cap, '-', 'Color', 'r');
plot(data.time, data.body, '-', 'Color', 'k');
plot(data.time, data.base, '-', 'Color', 'bl');
plot(data.time(ventInd), data.cap(ventInd), 'o');
plot(data.time(TRInd), data.cap(TRInd), 'o');
xlabel('Time [min]');
ylabel('Temperature [˚C]');
legend boxoff
legend('Cap', 'Body', 'Base', 'FontSize', 12, 'Location', 'NorthEast');
xlim([0, 100]);

subplot(2, 1, 2);
hold on
yyaxis left
plot(data.time, data.cap, '-', 'Color', 'bl');
xlabel('Time [min]');
ylabel('Temperature [˚C]');
yyaxis right
plot(data.time, data.crate, '-', 'Color', 'r');
plot(data.time(ventInd), data.crate(ventInd), 'o');
plot(data.time(TRInd), data.crate(TRInd), 'o');
set(gca, 'YScale', 'log');
ylabel('Heating rate [˚C/min]');
xlim([data.time(tPeak) - 15, data.time(tPeak) + 15]);

%% Open data.
dataOut = [ventTime, ventTemp, ventError, TRTime, TRTemp, TRError];
open dataOut;