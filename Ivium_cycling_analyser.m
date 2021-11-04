
%% Mark Buckwell, Electrochemical Innovation Lab, University College London
% Ivium_cycling_analyser, v1.0 2021

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

% Read in and process Ivium idf cycling files. Will split out all data
% types and save into separate files as length of each data type can vary
% depending on sampling parameters in method.

%% Set cycling analysis parameters.
VLims = [2.5 , 4.2]; % Cycling voltage limits, assuming cycling is reasonably 'normal'.
VDev = 0.02; % Deviation in V to class as intentional change (i.e. above noise).
nominalCapacity = 4200; % Nominal cell capacity to determine fade.

%% Read in files.
addpath(cd); % Add Matlab folder to path.
[vFiles, vLoc] = uigetfile('*.idf', 'Multiselect', 'On'); % Choose file.
cellName = textscan(vLoc, '%s', 'Delimiter', '\'); % Get data location.
cellName = cellName{1}{end}; % Get cell name.
importOpts = delimitedTextImportOptions("NumVariables", 6); % Create import options.
importOpts.Delimiter = ["\t", " "]; % Define delimiters.
cd(vLoc) % Set working folder to data folder.
if ischar(vFiles) == 1 % If a single file is selected.
    vFiles = {vFiles}; % Convert it to a single cell
end
cyclingData = zeros(1, 3); % Array to store voltage data.
VtTest = 0; % Initialise voltage time counter.
temp1Data = zeros(1, 2); % Array to store temp1 data.
T1tTest = 0; % Initialise temperature1 time counter.
temp2Data = zeros(1, 2); % Array to store temp2 data.
T2tTest = 0; % Initialise temperature2 time counter.
chargeData = zeros(1, 2); % Array to store charge data.
QtTest = 0; % Initialise charge time counter.
for iFiles = 1 : numel(vFiles) % Iterate through files.
    disp(vFiles{iFiles}); % Show file being processed.
    fileData = readtable([vLoc, vFiles{iFiles}], importOpts); % Read in data.
    emptyCells = cellfun(@isempty, fileData{:, :}); % Find empty cells.
    for iReplace = 1 : 6 % Iterate through first 6 rows, assuming these are the only ones containing useful data.
        
    end
    [~, cellID] = max(~cellfun(@isempty,...  % Find test ID.
        strfind(fileData.Var1, 'Title=')));
    cellID = textscan(fileData.Var1{cellID}, '%s', 'Delimiter','=');
    disp(cellID{1}{2}); % Display test ID to confirm cell ID is consistent.
    %% Find indices of each data type and split out.
    [~, VStart] = max(~cellfun(@isempty,... % Find start of voltage data.
        strfind(fileData.Var1, 'primary_data')));
    VStart = VStart + 3;
    [~, TStart] = max(~cellfun(@isempty,... % Find start of temperature data.
        strfind(fileData.Var1, 'analog_data')));
    VEnd = TStart - 2; % Align to correct row.
    TStart = TStart + 4; % Align to correct row.
    [~, QStart] = max(~cellfun(@isempty,... % Find start of charge data.
        strfind(fileData.Var1, 'Q_data')));
    TEnd = QStart - 2; % Align to correct row.
    QStart = QStart + 4; % Align to correct row.
    %% Porcess voltage data.
    VTime_min = str2double(fileData.Var1(VStart : VEnd)) ./ 60 + VtTest; % Get voltage time data.
    VtTest = VTime_min(end) + 0.1; % Store new voltage time counter value.
    DVoltage = str2double(fileData.Var4(VStart : VEnd)); % Get discharge voltage data.
    DVoltage(isnan(DVoltage)) = 0; % Set NaN values to zero.
    CVoltage = str2double(fileData.Var5(VStart : VEnd)); % Get charge voltage data.
    CVoltage(isnan(CVoltage)) = 0; % Set NaN values to zero.
    Voltage_V = DVoltage + CVoltage; % Sum to combine voltage values.
    Current_A = str2double(fileData.Var3(VStart : VEnd)); % Get current data.
    cyclingData = [cyclingData; VTime_min, Voltage_V, Current_A]; % Append voltage data.
    %% Process temperature data.
    allTemp = fileData.Var3(TStart : TEnd); % Get all temperature data.
    TMids = find(cellfun(@isempty, allTemp)); % Find middle of temperature data.
    T1_C = str2double(allTemp(1 : TMids(1) - 1)) .* 200; % Get first temperature data.
    T1Time_min = str2double(fileData.Var1(TStart : TStart + TMids(1) - 2)) ./ 60 + T1tTest; % Get temp time data.
    T1tTest = T1Time_min(end) + 0.1; % Store new temp time counter value.
    temp1Data = [temp1Data; T1Time_min, T1_C]; % Append temp1 data.
    T2_C = str2double(allTemp(TMids(2) + 1 : end)) .* 200; % Get second temperature data.
    T2Time_min = str2double(fileData.Var1(TStart + TMids(2) : TEnd)) ./ 60 + T2tTest; % Get temp time data.
    T2tTest = T2Time_min(end) + 0.1; % Store new temp time counter value.
    temp2Data = [temp2Data; T2Time_min, T2_C]; % Append temp1 data.
    %% Process charge data.
    QTime_min = str2double(fileData.Var1(QStart : end - 3)) ./ 60 + QtTest; % Get charge time data.
    QtTest = QTime_min(end) + 0.1; % Store new charge time counter value.
    Q_C = str2double(fileData.Var3(QStart : end - 3)); % Get charge data.
    chargeData = [chargeData; QTime_min, Q_C]; % Append charge data.
end
cyclingData(1, :) = []; % Remove starting voltage row.
temp1Data(1, :) = []; % Remove starting temp1 row.
temp2Data(1, :) = []; % Remove starting temp2 row.
chargeData(1, :) = []; % Remo starting charge row.
save(strcat(cellName, '_cycling.txt'), 'cyclingData', '-ascii'); % Save cycling data.
save(strcat(cellName, '_temp1.txt'), 'temp1Data', '-ascii'); % Save temp1 data.
save(strcat(cellName, '_temp2.txt'), 'temp2Data', '-ascii'); % Save temp2 data.
save(strcat(cellName, '_charge.txt'), 'chargeData', '-ascii'); % Save charge data.

%% Extract cycling from peaks in voltage.
VCharge = cyclingData(:, 2); % Get voltage data for charge.
VCharge(VCharge < VLims(2) - VDev) = VLims(1); % Remove artefacts above max cycle V.
VCharge(VCharge > VLims(2) + VDev) = VLims(1); % Remove data below max cycle V.
VCharge = -1 .* (VCharge - max(VCharge)); % Invert data.
[~, cLocs] = findpeaks(VCharge, 'MinPeakHeight', diff(VLims) - VDev);
cLocs = cLocs - 1; % Align timepoints to charge rather than discharge start/end.
VDischarge = cyclingData(:, 2); % Get voltage data for discharge.
VDischarge(VDischarge < VLims(1) - VDev) = VLims(2); % Remove artefacts below min cycle V.
VDischarge(VDischarge > VLims(1) + VDev) = VLims(2); % Remove data above min cycle V.
[~, dLocs] = findpeaks(VDischarge, 'MinPeakHeight', VLims(2) - VDev);
dLocs = dLocs - 1; % Align  timepoints to discharge rather than charge start/end.
dLocs(dLocs < cLocs(1)) = []; % Remove discharging prior to first charge.
for iCharge = 1 : numel(cLocs) - 1 % Remove erroneous charge timepoints.
    cDiff = max(cyclingData(cLocs(iCharge) : cLocs(iCharge + 1), 2)) -...
        min(cyclingData((iCharge) : cLocs(iCharge + 1), 2));
    if abs(diff(VLims) - cDiff) > (2 * VDev) % If large voltage range.
        cLocs(iCharge) = NaN; % Mark erroneous timepoint.
    end
end
for iDischarge = 1 : numel(dLocs) - 1 % Remove erroneous discharge timepoints.
    dDiff = max(cyclingData(dLocs(iDischarge) : dLocs(iDischarge + 1), 2)) -...
        min(cyclingData(dLocs(iDischarge) : dLocs(iDischarge + 1), 2));
    if abs(diff(VLims) - dDiff) > (2 * VDev) % If large voltage range.
        dLocs(iDischarge) = NaN; % Mark erroneous timepoint.
    end
end
cLocs(isnan(cLocs)) = []; % Remove erroneous timepoints.
dLocs(isnan(dLocs)) = []; % Remove erroneous timepoints.
if numel(cLocs) ~= numel(dLocs) % Alert  if unmatched cycles.
    disp(strcat('Unmatched cycling;', {' '},...
        string(numel(cLocs)), ' charging peaks and', {' '},...
        string(numel(dLocs)), ' discharging peaks found'));
    if dLocs(1) < cLocs(1) % This shouldn't happen.
        disp('First discharge prior to first charge.');
    end
else
    disp(strcat(string(numel(cLocs)), ' cycles detected.'));
end
capacityCharge = cellData{cLocs, 'Capacity_mAh'}; % Record charge capacities.
chargeTimes = cyclingData(cLocs, 1); % Record charge timepoints.
capacityDischarge = cellData{dLocs, 'Capacity_mAh'}; % Record discharge capacities.
dischargeTimes = cyclingData(dLocs, 1); % Record discharge timepoints.
% Assign cycle numbers based on charging and discharging indices.
for iCycle = 1 : numel(cLocs) - 1
    cellData{cLocs(iCycle) : cLocs(iCycle + 1) - 1,...
        'Cyc_num'} = iCycle; % Define first charge begins first cycle.
end
cellData{cLocs(end) : end, 'Cyc_num'} = iCycle + 1;

%% Plot capacity data.
figure;
subplot(2, 1, 1);
hold on
plot(cellData{:, 'Time_days'}, cellData{:, 'Capacity_mAh'});
xlabel('Time [days]');
ylabel('Capacity [mAh]');
subplot(2, 1, 2);
hold on
plot((capacityCharge ./ nominalCapacity) .* 100, 'o');
plot((capacityDischarge ./ nominalCapacity) .* 100, 'o');
ylim([80, 100]);
xlabel('Cycle number');
ylabel('Capacity [%]');
legend('Charge capacity', 'Discharge capacity', 'Location', 'NorthEast');
legend boxoff