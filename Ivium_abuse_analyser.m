function [voltageData, eventPoints] = Ivium_abuse_analyser()
% Read in and process Ivium idf files.
%   Will sort out the time, voltage, current and temperature data, as well
%   as find the event points, e.g. voltage drop(s) and peak temperature.
%   This will not process data to determine cycling, however.

set(groot, 'defaultFigureWindowStyle', 'Docked','defaultFigureColor', 'White',...
    'defaultAxesFontSize', 16, 'defaultAxesFontName', 'Arial',...
    'defaultLineLineWidth', 1.3, 'defaultAxesTickDir', 'Out',...
    'defaultAxesTickDirMode', 'Manual', 'defaultAxesLineWidth', 2,...
    'defaultFigureColormap', gray(64));

%% Read in file.
[vFile, vLoc] = uigetfile('*.idf', 'Multiselect', 'Off'); % Choose file.
importOpts = delimitedTextImportOptions("NumVariables", 6); % Create import options.
importOpts.Delimiter = ["\t", " "]; % Define delimiters.
allData = readtable([vLoc, vFile],importOpts); % Read in data.

%% Find indices of each data type and split out.
[~, VStart] = max(~cellfun(@isempty,... % Find start of voltage data.
    strfind(allData.Var1, 'primary_data')));
VStart = VStart + 3;
[~, TStart] = max(~cellfun(@isempty,... % Find start of temperature data.
    strfind(allData.Var1, 'analog_data')));
VEnd = TStart - 2; % Align to correct row.
TStart = TStart + 4; % Align to correct row.
[~, QStart] = max(~cellfun(@isempty,... % Find start of charge data.
    strfind(allData.Var1, 'Q_data')));
TEnd = QStart - 2; % Align to correct row.
QStart = QStart + 4; % Align to correct row.
Time_min = str2double(allData.Var1(VStart : VEnd)) ./ 60; % Get time data.
Voltage_V = str2double(allData.Var5(VStart : VEnd)); % Get voltage data.
Current_A = str2double(allData.Var3(VStart : VEnd)); % Get current data.
Q_C = str2double(allData.Var3(QStart : end - 3)); % Get charge data.
allTemp = allData.Var3(TStart : TEnd); % Get all temperature data.
TMids = find(cellfun(@isempty, allTemp)); % Find middle of temperature data.
T1_C = str2double(allTemp(1 : TMids(1) - 1)) .* 200; % Get first temperature data.
T2_C = str2double(allTemp(TMids(2) + 1 : end)) .* 200; % Get second temperature data.

%% Ask user for offsets, combine all data into a table and find change points.
disp(vFile); % Show file being processed.
tOffsets = inputdlg({'Voltage offset [s]:', 'ARC offset [s]:'},...
    'Offsets:', [1 35], {'0', '0'}); % Ask user for test time offsets.
Time_min = Time_min + double(tOffsets{1}) - double(tOffsets{2}); % Assume ARC starts at 0 min.
voltageData = table(Time_min, Voltage_V, Current_A, Q_C, T1_C, T2_C);
figure;
findpeaks(-1 .* (gradient(Voltage_V)), 'MinPeakHeight', 0.4); % Find drop peaks.
ylabel('|dV| [V]');
box off
grid off
hold on
yyaxis right
plot(Voltage_V);
ylabel('Voltage [V]');
[~, dropLocs] = findpeaks(-1 .* (gradient(Voltage_V)),...
    'NPeaks', 2, 'MinPeakHeight', 0.4); % Find up to 2 peaks.
dropLocs = dropLocs - 1; % Align to timepoints of drop.
eventPoints = Time_min(dropLocs); % Output drop times.
end