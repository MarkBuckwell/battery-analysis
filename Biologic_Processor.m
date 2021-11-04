
%% Mark Buckwell, Electrochemical Innovation Lab, University College London
% Biologic_Cycling_Processor, v1.0 2021

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

%% This script will take a messy set of Biologic cycling files (messy being
% inconsistent measurement setups, lengths, test conditions etc) and find
% the charge/discharge cycles in each file, putting them together into a
% full set of cycle data. It will also extract the EIS data and sort it
% into appropriate locations in the cycle life of the cell, as well as
% perform a consistent analysis on each measurement separately. Note that
% Zfit is required for EIS analysis.

% The script will ignore files that don't contain data or are invalid, but
% will notify the user which files these are. So, the ideal
% input is a folder of files corresponding to the cycle life of a single
% cell.

% The script behaves by joining together all the data recorded for a given
% cell/sample, then sorting it by the time point for each data point. The
% cycles are then determined by assessing peaks in the normalised voltage
% divided by the normalised current.

% The time format in the file can be the test time or real time, both will
% be processed and used to sort the full dataset. However, it is assumed
% that the same time format is used for all files.

%% IMPORTANT INFORMATION ON DATA SORTING.
% Sorting of data in the correct time order will only work if the data
% files are labelled correctly, wherein the first potion of the name is the
% date, in the format YYMMDDHHMM, or YYMMDD (which will get padded to
% 00:00), followed by an underscore. E.g. all files for a single test will
% have the correct start time, and each subsequent set of test files will
% be shifted to after the previous file. However, this will also work if,
% for examples, the file set are numbered 1_, 2_, 3_ etc, though they will
% still be padded.
%% READ ABOVE.

%% Initialise parameters and counters.
VLims = [2.5 , 4.2]; % Cycling voltage limits, assuming cycling is reasonably 'normal'.
VDev = 0.02; % Deviation in V to class as intentional change (i.e. above noise).
nominalCapacity = 4200; % Nominal cell capacity to determine fade.
VEIS = [0.015, 4.1]; % Cutoff voltage for [GEIS, PEIS] in mV.
IEIS= [140, 300]; % Cutoff current for [GEIS, PEIS] in mA.
cellHeaders = {'Time_days', 'Voltage_V', 'Temp_C', 'Current_A',...
    'Capacity_mAh', 'Freq_Hz', 'Z_Re', 'Z_Im', 'Cyc_num'};
% More may be added as required, indexing and processing of eisting
% variables should still work in new variables are added in AFTER existing
% ones. Otherwise, need to further generalise this script for importing
% columns of data.

%% Choose folder of files to work through.
cd(uigetdir);
fileDir = dir;
warning('Off'); % Ignore warnings.
try % To read cell capacity from setting file.
    [~, ~, fExt] = fileparts(fileDir(3).name); % Get file extension.
    fileChoice = fopen(fileDir(3).name); % Open file.
    cellChar = textscan(fileChoice, '%s', 'Delimiter', '\n');
    cellCap = find(contains(cellChar{1}, 'Battery capacity'));
    cellCap = textscan(cellChar{1}{cellCap}, '%s', 'Delimiter', ':');
    cellCap = textscan(cellCap{1}{2}, '%f', 'Delimiter', ' ');
    cellCap = cellCap{1} .* 1000; % User input capacity in mAh.
    if cellCap > 0 % Check that value exists.
        nominalCapacity = cellCap; % Overwrite value in script.
    end
catch
    disp('Invalid capacity value given.');
    disp(' ');
end

%% Iterate through all files and collate all data into single table.
nVariables = numel(cellHeaders);
cellData = array2table(zeros(1, nVariables), 'VariableNames', cellHeaders);
currentDate = 0; % Initialise the file date counter.
realTime = 0; % Initialise time format as test time.
realTimes = {'0'}; % Initialise real time counter.
for iFiles = 1 : numel(fileDir)
    [~, ~, fExt] = fileparts(fileDir(iFiles).name); % Get file extension.
    if fileDir(iFiles).isdir == 0 && ~contains('.txt', fExt) == 0
        disp(fileDir(iFiles).name); % State file name.
        try % Import data to process.
            fileChoice = fopen(fileDir(iFiles).name); % Open file.
            headers = textscan(fgetl(fileChoice), '%s', 'Delimiter', '\t');
            if contains(headers{1}, 'BT-Lab ASCII FILE') == 1 % For headed files.
                headSkip = textscan(fgetl(fileChoice),...
                    '%s', 'Delimiter', ':');
                headSkip = str2double(headSkip{1}{end}); % Get header lines.
                for iSkip = 1 : headSkip - 3 % Skip through to data headers.
                    fgetl(fileChoice);
                end
                headers = textscan(fgetl(fileChoice), '%s', 'Delimiter', '\t');
            end
            headers = erase(headers{1},... % Sort out headers.
                ['/', {' '}, '.', '-', '(', ')', '|', '°', 'µ', '%']);
            dataTable = array2table(textscan(fileChoice,...
                repmat('%s', 1, numel(headers)), 'Delimiter', '\t'),...
                'VariableNames', headers); % Get data.
            for iConvert = 1 : numel(headers) % Convert to proper formats.
                if contains(dataTable{1, iConvert}{1}(1), '/') == 1
                    realTime = 1; % Set time format to realtime and record separately.
                    realTimes = [realTimes; dataTable{1, iConvert}{1}];
                    dataTable{1, iConvert}{1} = zeros(numel(...
                        dataTable{1, iConvert}{1}), 1); % Replace with 0s.
                else % Convert cells to double.
                    dataTable{:, iConvert}{1} =...
                        str2double(dataTable{:, iConvert}{1});
                end
            end
            
            %% Check file ordering if using test time format.
            % There could be an issue if multiple files are made on the
            % same date but have erroneous date labels.
            if realTime == 0 % If date format is test time
                fileDate = textscan(fileDir(iFiles).name, '%d', 'Delimiter', '_');
                dateLength = numel(num2str(fileDate{1})); % Get length of date.
                fileDate = int32(fileDate{1} .* (10 ^ (10 - dateLength))); % Add trailing zeros.
                if fileDate > currentDate % Adjust file timepoints.
                    timeOffset = cellData{end, 'Time_days'} + 10; % Add 10 second gap.
                    dataTable.times{1} = dataTable.times{1} + timeOffset;
                    currentDate = fileDate; % Store new date.
                elseif fileDate == currentDate % Use existing shift for tests on same date.
                    dataTable.times{1} = dataTable.times{1} + timeOffset;
                else
                    disp('Date calibration error.');
                end
            end
            
            %% Concatenate data.
            try % Case that all data types are present.
                cellData = [cellData ; array2table([dataTable.times{1},...
                    dataTable.EcellV{1}, dataTable.TemperatureC{1},...
                    dataTable.ImA{1}, dataTable.CapacitymAh{1},...
                    dataTable.freqHz{1}, dataTable.ReZOhm{1},...
                    dataTable.ImZOhm{1}, zeros(fileLength, 1)],...
                    'VariableNames', cellHeaders)];
                disp('File contains all data types.');
            catch
                fileLength = size(dataTable.EcellV{1}, 1);
                try % Case that no EIS data is present.
                cellData = [cellData ; array2table([dataTable.times{1},...
                    dataTable.EcellV{1}, dataTable.TemperatureC{1},...
                    dataTable.ImA{1}, dataTable.CapacitymAh{1},...
                    nan(fileLength, 3), zeros(fileLength, 1)],...
                    'VariableNames', cellHeaders)];
                disp('File contains only cycling data.');
                catch
                    try % Case that only EIS data is present.
                    cellData = [cellData ; array2table([dataTable.times{1},...
                        dataTable.EcellV{1}, nan(fileLength, 1),...
                        dataTable.ImA{1}, nan(fileLength, 1),...
                        dataTable.freqHz{1}, dataTable.ReZOhm{1},...
                        dataTable.ImZOhm{1}, zeros(fileLength, 1)],...
                        'VariableNames', cellHeaders)];
                    disp('File contains only EIS data.');
                    catch
                        try % Case that only resting data is present.
                        cellData = [cellData ; array2table([dataTable.times{1},...
                            dataTable.EcellV{1}, nan(fileLength, 6),...
                            zeros(fileLength, 1)], 'VariableNames', cellHeaders)];
                        disp('File contains only rest data.');
                        catch
                            ('Data read error.');
                        end
                    end
                end
            end
        catch
            disp('File read error.');
        end
        disp(' ');
    end
end
disp('All files now checked.');
cellData(1, :) = []; % Remove first row of zeros.
realTimes(1) = []; % Remove initial zero time.
if realTime == 1 % Replace and convert to datetime.
    cellData.Time_days = datetime(realTimes, 'Format', 'preserveinput' );
else % Convert test time to days.
    cellData{:, 'Time_days'} = cellData{:, 'Time_days'} ./ (60 * 60 * 24);
end
cellData = sortrows(cellData); % Put all data in chronological order.
cellData{:, 'Capacity_mAh'} = fillmissing(... % Replace NaN capacity values.
    cellData{:, 'Capacity_mAh'}, 'Nearest');
cellData{:, 'Current_A'} = cellData{:, 'Current_A'} ./ 1000; % Convert current to A.

%% Extract cycling from peaks in voltage.
VCharge = cellData{:, 'Voltage_V'}; % Get voltage data for charge
VCharge(VCharge < VLims(2) - VDev) = VLims(1); % Remove artefacts above max cycle V.
VCharge(VCharge > VLims(2) + VDev) = VLims(1); % Remove data below max cycle V.
VCharge = -1 .* (VCharge - max(VCharge)); % Invert data.
[~, cLocs] = findpeaks(VCharge, 'MinPeakHeight', diff(VLims) - VDev);
cLocs = cLocs - 1; % Align timepoints to charge rather than discharge start/end.
VDischarge = cellData{:, 'Voltage_V'}; % Get voltage data for discharge.
VDischarge(VDischarge < VLims(1) - VDev) = VLims(2); % Remove artefacts below min cycle V.
VDischarge(VDischarge > VLims(1) + VDev) = VLims(2); % Remove data above min cycle V.
[~, dLocs] = findpeaks(VDischarge, 'MinPeakHeight', VLims(2) - VDev);
dLocs = dLocs - 1; % Align  timepoints to discharge rather than charge start/end.
dLocs(dLocs < cLocs(1)) = []; % Remove discharging prior to first charge.
for iCharge = 1 : numel(cLocs) - 1 % Remove erroneous charge timepoints.
    cDiff = max(cellData{cLocs(iCharge) : cLocs(iCharge + 1), 'Voltage_V'}) -...
        min(cellData{cLocs(iCharge) : cLocs(iCharge + 1), 'Voltage_V'});
    if abs(diff(VLims) - cDiff) > (2 * VDev) % If large voltage range.
        cLocs(iCharge) = NaN; % Mark erroneous timepoint.
    end
end
for iDischarge = 1 : numel(dLocs) - 1 % Remove erroneous discharge timepoints.
    dDiff = max(cellData{dLocs(iDischarge) : dLocs(iDischarge + 1), 'Voltage_V'}) -...
        min(cellData{dLocs(iDischarge) : dLocs(iDischarge + 1), 'Voltage_V'});
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
chargeTimes = cellData{cLocs, 'Time_days'}; % Record charge timepoints.
capacityDischarge = cellData{dLocs, 'Capacity_mAh'}; % Record discharge capacities.
dischargeTimes = cellData{dLocs, 'Time_days'}; % Record discharge timepoints.
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

%% Plot temperature data.
cellT = cellData{:, 'Temp_C'};
cellT(cellT > 100) = mode(cellT); % Remove errors due to thermocouple unplugging.
figure;
% 30/12 good for 1C.
% 50/50 good for >= 3C.
try
    plot(cellData{:, 'Time_days'}, cellT);
    xlabel('Time [days]');
    ylabel(['Temperature [' char(176) 'C]']);
    box off
    I = cellData{:, 'Current_A'};
    I(isnan(I)) = 0;
    V = cellData{:, 'Voltage_V'};
    V(V > 4.25) = 4.2;
    V(V < 2.45) = 2.5;
    V(isnan(V)) = 0;
    t = cellData{:, 'Time_days'} .* 24 .* 60 .* 60;
    thermWork = 0.066 .* (cellT - 22.5) * 1619;
    elecP = abs(I) .* V;
    elecWork = gradient(t) .* elecP;
    workLoss = elecWork - thermWork;
catch
    disp('No cycling available for plotting temperature.');
end

%% Extract EIS measurements using peaks in the frequency data.
%  Also corrected for large jumps i.e. if cell setup has been changed due
%  to detachment/reattachment.
freqData = cellData{:, 'Freq_Hz'}; % Get frequency data.
freqData(isnan(freqData)) = 0; % Replace NaN values.
% Find EIS measurements using peaks for onset of test frequency range.
[~, ZStartInds] = findpeaks(freqData, 'MinPeakProminence', 8000);
[~, ZEndInds] = findpeaks(-1 .* freqData, 'MinPeakProminence', 8000);
ZEndInds = ZEndInds - 1; % Correct indices to end of frequency sweep.
ZEndInds(ZEndInds < ZStartInds(1)) = []; % Remove any end indices before first start.
if numel(ZEndInds) < numel(ZStartInds) % If the final test has no end peak.
    ZEndInds = [ZEndInds; numel(freqData)]; % Use the final data index.
end
impedanceData = zeros(numel(ZStartInds), 3); % Matrix to store EIS metrics.
impedanceTests = repmat({''}, numel(ZStartInds), 5); % Array to store measurements.
for iEIS = 1 : numel(ZStartInds) % For each individual test.
    impedanceTests{iEIS, 1} = cellData{ZStartInds(iEIS) : ZEndInds(iEIS), 'Freq_Hz'};
    impedanceTests{iEIS, 2} = cellData{ZStartInds(iEIS) : ZEndInds(iEIS), 'Z_Re'};
    impedanceTests{iEIS, 3} = cellData{ZStartInds(iEIS) : ZEndInds(iEIS), 'Z_Im'};
    impedanceTests{iEIS, 4} = cellData{ZStartInds(iEIS) : ZEndInds(iEIS), 'Voltage_V'};
    impedanceTests{iEIS, 5} = cellData{ZStartInds(iEIS) : ZEndInds(iEIS), 'Current_A'};
    impedanceData(iEIS, 1) = cellData{ZStartInds(iEIS),...
        'Cyc_num'}; % Get cycle number.
    [~, iRzc] = min(abs(impedanceTests{iEIS, 3})); % Get zero-crossing resistance.
    impedanceData(iEIS, 2) = impedanceTests{iEIS, 2}(iRzc);
    fzc(iEIS) = impedanceTests{iEIS, 1}(iRzc);
end

% Determine presence of disconnection jumps using zero-corssing resistance.
% Assumes that the difference between consecutive points is fairly
% consistent, so a jump may be corrected by removing the difference between
% points and adding the mean difference between all previous points.
figure;
subplot(2, 1, 2); % Plot zero crossing resistance.
plot(impedanceData(:, 1), impedanceData(:, 2), 'o');
xlabel('Cycle number');
ylabel('R_{zc} [\Omega]');
box off
hold on
topPlot = subplot(2, 1, 1);
disHeight = 0.0006; % Starting peak height
findpeaks(abs(...
    gradient(impedanceData(:, 2)) ./...
    gradient(impedanceData(:, 1))),...
    'MinPeakHeight', disHeight); % Trial peaks of R zc spikes for disconnection.
disHeight = str2double(inputdlg({'Enter threshold height:'},...
    'Disconnection inspection', [1, 35], {char(string(disHeight))}));
[~, disLocs] = findpeaks(abs(...
    gradient(impedanceData(:, 2)) ./...
    gradient(impedanceData(:, 1))), 'MinPeakWidth', 1,...
    'MinPeakHeight', disHeight); % Get peaks of R zc spikes for disconnection.
delete(topPlot);
for iDis = 1 : numel(disLocs) % Start with mean and std of difference between vals.
    meanZCDiff = mean(diff(impedanceData(1 : disLocs(iDis) - 1, 2)));
    stdZCDiff = std(diff(impedanceData(1 : disLocs(iDis) - 1, 2)));
    for iCorr = disLocs(iDis) : numel(ZStartInds) % Correct tests.
        impedanceTests{iCorr, 2} = impedanceTests{iCorr, 2} -...
            (impedanceData(disLocs(iDis), 2) -...
            impedanceData(disLocs(iDis) - 1, 2)) + meanZCDiff;
    end
    impedanceData(disLocs(iDis) : end, 2) =... % Correct using mean difference.
    impedanceData(disLocs(iDis) : end , 2) -...
    (impedanceData(disLocs(iDis), 2) -...
    impedanceData(disLocs(iDis) - 1, 2)) + meanZCDiff;
end

subplot(2, 1, 2); % Plot, fit and calibrate EIS data.
plot(impedanceData(:, 1), impedanceData(:, 2), 'x');
zcType = fittype('a + (b * (log(c * x)))',... % Logarithmic fit
    'dependent', {'y'}, 'independent', {'x'},...
    'coefficients',{'a', 'b', 'c'});
zcFit = fit(impedanceData(2 : end, 1), impedanceData(2 : end, 2), 'poly1');
% zcFit = fit(impedanceData(2 : end, 1), impedanceData(2 : end, 2), zcType,...
%     'Lower', [0, 0, 0], 'Upper', [0.1, 0.1, 5]);
plot(impedanceData(2 : end, 1), zcFit(impedanceData(2 : end, 1)),...
    '-', 'Color', 'k');
legend('Raw data', 'Corrected', 'Fit', 'Location', 'SouthEast');
legend boxoff
subplot(2, 1, 1); % Plot corrected Bodes.
xlabel('Z_{Re} [\Omega]');
ylabel('Z_{Im} [\Omega]');
hold on
for iEIS = 1 : numel(ZStartInds) % For each individual test.
    plot(impedanceTests{iEIS, 2}, impedanceTests{iEIS, 3});
end

%% Batch analyse EIS data using Zfit.
figure;
EISFits = repmat({''}, size(impedanceTests, 1), 2); % To receive fitting params.
trialCircuit = 's(L1,R1,p(R1,E2),p(s(R1,H2),E2))'; % Circuit to try fitting.
testParams =  [2.72e-07,-4.21e-03,8.73e-03,1.41e+02,-1.06e+00,2.89e-04,-1.75e-01,2.56e-01,2.92e+01,2.24e-01]; % Trial fitting parameters.
fitUpper = [Inf, Inf, Inf, Inf,Inf, Inf, Inf,Inf Inf, Inf]; % Fitting upper bounds.
fitLower = [0, 0, 0, -Inf,-Inf, 0.0001 -Inf,-Inf, -Inf, -Inf]; % Fitting lower bounds.
for iEIS = 1 : size(impedanceTests, 1) % Iterate through measurements.
    [~, fitStart] = min(abs(impedanceTests{iEIS, 3})); % Fit from RZC.
    fitEnd = numel(impedanceTests{iEIS, 3}) - 10; % Fit to 10th from end.
    [bestParams, ZBest] = Zfit([impedanceTests{iEIS, 1},...
        impedanceTests{iEIS, 2}, -1 .* impedanceTests{iEIS, 3}], 'z', ...
        trialCircuit, testParams, [fitStart : fitEnd], 'fitNP',...
        fitLower, fitUpper); % Try fitting data.
    [bestParams, ZBest] = Zfit([impedanceTests{iEIS, 1},...
        impedanceTests{iEIS, 2}, -1 .* impedanceTests{iEIS, 3}], 'z', ...
        trialCircuit, bestParams, [fitStart : fitEnd], 'fitNP',...
        fitLower, fitUpper); % Use best params again and store new output.
    EISFits{iEIS, 1} = bestParams;
    EISFits{iEIS, 2} = ZBest;
end

%% Plot ZFit parameters.
figArrange = ceil(sqrt(numel(testParams))); % To arrange plotting.
varPos = find(isstrprop(trialCircuit,'Upper')); % Get variable positions.
varLabels = repmat({''}, 1, numel(testParams)); % Array for variable names.
iSkip = 0; % Initialse cell skipping index for double variable names.
for iLabels = 1 : numel(varPos)
    varLabels{iLabels + iSkip} = trialCircuit(varPos(iLabels) : varPos(iLabels) + 1);
    if varLabels{iLabels + iSkip}(2) == '2' % If double variable name.
        varLabels{iLabels + iSkip + 1} = varLabels{iLabels + iSkip};
        iSkip = iSkip + 1; % Skip and duplicate label.
    end
end
figure;
for iParams = 1 : numel(testParams)
    subplot(figArrange, figArrange, iParams);
    title(strcat(string(iParams), ':', varLabels{iParams}));
    hold on
    box off
    for iFits = 1 : size(EISFits, 1)
        plot(impedanceData(iFits, 1), EISFits{iFits, 1}(iParams), 'o',...
            'Color', 'r');
    end
end
    

%% Determine break points when cell removed for scanning. IN PROGRESS.
noCur = cellData{:, 'Current_A'}; % Get current data.
noCur(isnan(noCur)) = 0; % Set NaN current to 0.
restV = cellData{:, 'Voltage_V'}; % Get voltage data.
restV(abs(noCur) > 0.01) = 0; % Set rest periods to 0.

%% Explore other plotting options.
figure;
subplot(2, 1, 1);
plot(nominalCapacity - capacityCharge(impedanceData(2 : end, 1)),...
    (impedanceData(2 : end, 2) - impedanceData(2, 2)) .* 1000, 'o');
hold on
box off
xlabel('Capacity loss [mAh]');
ylabel('\DeltaR_{zc} [m\Omega]');
xlim([0, 1000]);
subplot(2, 1, 2);
plot(impedanceData(2, 2) .* (nominalCapacity -...
    capacityCharge(impedanceData(2 : end, 1))),...
    impedanceData(2 : end, 2), 'o');
hold on
box off
xlabel('R0(Cnom-Ci)');
ylabel('Ri');