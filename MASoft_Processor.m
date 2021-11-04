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

%% Import mass spectroscopy file.
% Set up to accept Hiden RGA MS .csv files, wherein each spectrum is a row
% and the rows are sequential spectra taken over time.
cal_data = readtable('mass_spec_cal.xlsx'); % Read in calibration file.
[fName, fLoc] = uigetfile('*.csv'); % Choose file to analyse.
cd(fLoc);
fData = readmatrix(strcat(fLoc, fName), 'NumHeaderLines', 22);
if sum(isnan(fData(:, end))) > 0 % If final spectrum is incompleted.
    fData(end, :) = []; % Delete final spectrum.
end
masses = fData(1, 4 : end); % Get mass axis values.
times = fData(2 : end, 3) ./ (1000 * 60 * 60); % Get time axis values.
fData(:, 1 : 3) = []; % Remove unused data.
fData(1, : ) = []; % Remove unused data.

%% Initial cropping of data to reduce memory usage.
figure;
semilogy(fData); % Show all data.
xlabel('Sample index');
ylabel('Counts [a.u.]');
box off
cropParams = {'1', num2str(size(fData, 1))};
cropParams = inputdlg({'Start index:', 'End index:'},...
    'Crop parameters', [1 35], cropParams); % User refinement.
fData = fData(str2double(cropParams{1}) : str2double(cropParams{2}), :); % Crop data.
times = times(str2double(cropParams{1}) : str2double(cropParams{2}), :); % Crop time.

%% Plot summed spectrum for whole measurement and find masses present.
sumSpec = sum(fData, 1); % Sum data in time for each mass value.
sumSpec(sumSpec <= 0) = min(abs(sumSpec)); % Remove zeroes and negative values.
subplot(2, 2, 1);
semilogy(masses, fData(1, :)); % Show a spectrum.
hold on
plot(masses, (0 .* masses) + 1E-12);
xlabel('Mass [amu]');
ylabel('Intensity [a.u.]');
box off
ylim([min(fData(1, :)), max(fData(1, :))]);
subplot(2, 2, 2); % Determine noise level.
noisePlot = plot(masses, fData(1,:));
xlabel('Mass [amu]');
ylabel('Intensity [a.u.]');
ylim([1E-11, 5E-10]);
box off
inspectPeaks = {'3e-8', '3E-9', '0.1', '0.5', '4E-11'}; % Initial inspection parameters.
% Good levels: 8.3e-11, 1.4E-10
iFind = 1; % Initialise peak inspection index.
warning off
while iFind == 1 % Allow user to adjust peak finding parameters.
    hold on
    noiseThresh2 = plot(masses, (0 .* masses) + str2double(inspectPeaks{5}));
    peakParams = str2double(inspectPeaks);
    subplot(2, 2, 3 : 4); % Identify mass peaks.
    hold off
    findpeaks(sumSpec, masses, 'MinPeakHeight', peakParams(1),...
        'MinPeakProminence', peakParams(2),...
        'MinPeakWidth', peakParams(3),...
        'MinPeakDistance', peakParams(4)); % Look for peaks.
    set(gca, 'YScale', 'Log');
    xlabel('Mass [amu]');
    ylabel('Intensity [a.u.]');
    ylim([1E-9, 1E-4]);
    box off
    acceptPeaks = questdlg('Accept peak finding?', 'Options',...
        'Accept', 'Adjust parameters', 'End', 'Adjust parameters');
    switch acceptPeaks
        case 'Accept'
            [mPks, mLocs, mWidths] = findpeaks(sumSpec, masses,...
                'MinPeakHeight', peakParams(1),...
                'MinPeakProminence', peakParams(2),...
                'MinPeakWidth', peakParams(3),...
                'MinPeakDistance', peakParams(4)); % Commit peak finding.
            iFind = 0; % Reset peak finding index.
        case 'Adjust parameters' % Allow user to revise peak finding.
            inspectPeaks = inputdlg({'Min peak height:', 'Min peak prominence:',...
                'Min peak width;', 'Min peak distance:', 'Noise floor:'},...
                'Peak parameters', [1 35], inspectPeaks); % User refinement.
        case 'End'
            return % End the script.
    end
    subplot(2, 2, 2);
    delete(noiseThresh2); % Remove noise threshold marker.
end

%% Use identified peaks to track each species over time.
trackData = zeros(size(fData, 1), numel(mLocs)); % Matrix to store mass data.
clearData = fData; % Create new matrix for data processing.
clearData(clearData < str2double(inspectPeaks{5})) = 0; % Remove noise.
delete(noisePlot); % Remove previous plot.
clear showPeaks % Remove data if re-running this section.
clear spreadPeaks % Remove data if re-running this section
hold off
sampleSpread = 3; % Number of samples to include in each direction in spread data.
for iMass = 1 : numel(mLocs) % Iterate through masses identified.
    iShow = 1; % Counter to show some peaks.
	[~, mLower] = min(abs(masses -... % Trial lower bound using found width.
        (mLocs(iMass) - (mWidths(iMass) / 1) - 0)));
    [~, mUpper] = min(abs(masses -... % Trial upper bound using found width.
        (mLocs(iMass) + (mWidths(iMass) / 1) + 0)));
    for iSample = 1 : size(fData, 1) % Iterate through spectra in time.
        peakCrop = clearData(iSample, mLower : mUpper); % Crop out mass peak.
        massCrop = masses(mLower : mUpper); % Crop out mass range.
        if mLower < sampleSpread + 1 % If checking mass peak near start of range.
            mSpreadL = 1; % Set lower bound to first mass value.
        else % If checking mass away from start of range.
            mSpreadL = mLower - sampleSpread; % Set broader lower value.
        end
        if mUpper > numel(masses) - sampleSpread % If checking mass peak near end of range.
            mSpreadU = numel(masses); % Set upper bound to final mass value.
        else % If checking mass away from end of range.
            mSpreadU = mUpper + sampleSpread; % Set broader upper bound.
        end
        peakSpread = clearData(iSample, mSpreadL : mSpreadU); % Capture broader spectral range.
        massSpread = masses(mSpreadL : mSpreadU); % Capture broader mass range.
        clearData(iSample, mLower : mUpper) = 0; % Remove used data.
        if rem(iSample, 1) == 0 % Record and plot peaks.
            showPeaks{iShow, iMass} = [massCrop', peakCrop'];
            spreadPeaks{iShow, iMass} = [massSpread', peakSpread'];
            iShow = iShow + 1; % Increment showing index,
%             plot(masses(mLower : mUpper), peakCrop);
%             xlabel('Mass [amu]');
%             ylabel('Intensity [a.u.]');
%             title(strcat(string(iSample), '/', string(size(fData, 1))));
%             box off
%             pause(0.0001);
        end
        if numel(massCrop) == 1 && numel(peakCrop) == 1
            trackData = mean(diff(masses)) * peakCrop; % Average for single point.
        elseif numel(massCrop) > 1 && numel(peakCrop) > 1
            trackData(iSample, iMass) = trapz(massCrop, peakCrop); % Integrate peak.
        else
            disp('Integration error.');
        end
    end
end
disp('Counted all mass peaks.');

%% Check validity of tracked peak data.
[~, NInd] = min(abs(mLocs - 28)); % Find N2 peak index.
[~, pInd] = max(gradient(trackData(:, NInd))); % Find thermal runaway peak.
figure;
subplot(2, 1, 2);
semilogy([times(pInd), times(pInd)], [1e-9, 3e-9], 'Color', 'k');
hold on
box off
xlabel('Time [hrs]');
ylabel('Counts [a.u.]');
iCheck = 1; % Initialise checking index.
peakColours = jet(size(showPeaks, 1)); % Generate set of colours.
while iCheck > 0
    subplot(2, 1, 1);
    cla
    hold on
    for iShow = 1 : size(showPeaks, 1) % Show integrated peaks.
        plot(spreadPeaks{iShow, iCheck}(:,1), spreadPeaks{iShow, iCheck}(:,2)',...
            '-', 'Color', peakColours(iShow, 1 : 3));
        plot(showPeaks{iShow, iCheck}(:,1), showPeaks{iShow, iCheck}(:,2)',...
            'o', 'Color', 'k');
    end
    title('Blue: first scan. Red: final scan. Markers: integrated data.');
    set(gca, 'YScale', 'log');
    xlabel('Mass [amu]');
    ylabel('Intensity [a.u.]');
    subplot(2, 1, 2);
    checkTrackO = semilogy(times, trackData(:, iCheck), 'o');
    checkTrackL = semilogy(times, trackData(:, iCheck), '-');
    legend(num2str(mLocs(iCheck)), 'Location', 'NorthEast');
    legend boxoff
    keepPeak = questdlg('Keep tracked peak?', 'Options',...
        'Keep', 'Remove', 'End', 'Keep');
    switch keepPeak
        case 'Keep' % Keep tracked peak in dataset.
            iCheck = iCheck + 1; % Increment checking index.
        case 'Remove' % Remove tracked peak from dataset.
            mLocs(iCheck) = []; % Remove peak mass.
            trackData(:, iCheck) = []; % Remove peak data.
            showPeaks(:, iCheck) = []; % Remove integrated peak data.
            spreadPeaks(:, iCheck) = []; % Remove broader peak region.
        case 'End' % Accept all data from this peak.
            disp('All subsequent data accepted.');
            iCheck = 0; % End checking.
    end
    if iCheck > numel(mLocs) % Case that all peaks have been checked.
        iCheck = 0; % End checking loop.
    end
    delete(checkTrackO); % Remove plotted data.
    delete(checkTrackL); % Remove plotted data.
end

%% Crop data to period of interest.
[~, NInd] = min(abs(mLocs - 28)); % Find N2 peak index.
figure;
subplot(2, 1, 1); % Plot full spectrum window.
semilogy(times, trackData(:, NInd));
xlabel('Time [hr]');
ylabel('Intensity [a.u.]');
box off
subplot(2, 1, 2); % Plot cropped spectrum window.
hold on
xlabel('Time [hr]');
ylabel('Intensity [a.u.]');
iCrop = 1; % Initialise cropping index.
defineCrop = {'5.8', '8'}; % Initial cropping window.
while iCrop == 1
    defineCrop = inputdlg({'Start time [hrs]:', 'End time[hrs]:'}, 'Crop window',...
        [1 35], defineCrop); % User refinement.
    cropWindow = str2double(defineCrop);
    [~, cropLow] = min(abs(times - str2double(defineCrop{1}))); % Lower bound.
    [~, cropUp] = min(abs(times - str2double(defineCrop{2}))); % Upper bound.
    try
        delete(cropPlot);
    catch
    end
    cropPlot = plot(times(cropLow : cropUp), trackData(cropLow : cropUp, NInd));
    acceptCrop = questdlg('Accept crop window?', 'Options',...
        'Accept', 'Adjust window', 'End', 'Adjust window');
    switch acceptCrop
        case 'Accept' % Commit cropping.
            cropData = trackData(cropLow : cropUp, :); % Crop data.
            cropTime = times(cropLow : cropUp); % Crop time.
            iCrop = 0; % Reset peak finding index.
        case 'End'
            return % End the script.
    end
end
dataOut = [cropTime, cropData]; % Combine data.
dataOut = [NaN, mLocs; dataOut]; % Add masses. 
save(strcat(fName(1 : end - 4), '_window.txt'), 'dataOut', '-ascii');

%% Label masses and determine fractional contributions.
mFound = round(mLocs); % Round found masses to integers.
% Define reference masses (expected or known).
mRefs = [1, 14, 16, 16, 16, 17, 18, 20,...
    28, 28, 28, 29, 32, 34, 40, 44];
% Define lables to assign to each mass. Can have multiple per mass.
mLabels = {'H', 'N', 'O_{2}', 'H_{2}O', 'CO', 'H_{2}O', 'H_{2}O', 'Ar'...
    'N_{2}', 'CO', 'CO_{2}', 'N_{2}', 'O_{2}', 'O_{2}', 'Ar', 'CO_{2}'};
% Define cracking pattern ratios for each mass assignment.
mCrack = [10, 7.2, 11.4, 1.1, 0.9, 23, 100, 10.7,...
    100, 100, 11.4, 0.8, 100, 0.4, 100, 100];
% Define relative sensitivities for each mass assignment.
mSens = [0.44, 1, 0.86, 0.9, 1.05, 0.9, 0.9, 1.2,...
    1, 1.05, 1.4, 1, 0.86, 0.86, 1.2, 1.4];
scaleData = trackData; % Matrix for storing scaled/processed data.

%% Carry out a trial background subtraction.
try
    try % To extract the argon trace.
        [~, ArInd] = min(abs(mFound - 40)); % Find Ar peak index.
        ArBack = trackData(:, ArInd);
    catch
        disp('No argon data found.');
    end
    try % To extract the nitrogen trace.
        [~, NInd] = min(abs(mFound - 28)); % Find N peak index.
        NBack = trackData(:, NInd); %
    catch
        disp('No nitrogen data found.');
    end
    normComp = (NBack ./ max(NBack)) ./ (ArBack ./ max(ArBack)); % Create comparison data.
%     ArFilt = lowpass(ArBack, 0.01); % Filter to get background.
    ArFilt = smooth(ArBack, 10, 'rlowess');
    ArFilt(1 : 6) = ArFilt(7); % Remove filtering artefact.
    ArFilt(end - 5 : end) = ArFilt(end - 6); % Remove filtering artefact.
    figure;
    subplot(2, 2, 1);
    semilogy(times, ArBack ./ max(ArBack));
    hold on
    semilogy(times, NBack ./ max(NBack));
    box off
    legend('Argon', 'Nitrogen', 'Location', 'SouthEast');
    legend box off
    xlabel('Time [hours]');
    ylabel('Norm. counts');
    subplot(2, 2, 2);
    semilogy(times, normComp);
    hold on
    box off;
    xlabel('Time [hours]');
    ylabel('Norm. ratio');
    subplot(2, 2, 3);
    semilogy(times, ArBack);
    hold on
    semilogy(times, ArFilt);
    box off;
    xlabel('Time [hours]');
    ylabel('Ar counts [a.u.]');
    legend('Raw', 'Filtered', 'Location', 'SouthEast');
    legend box off
    subplot(2, 2, 4);
    try % To extract test trace.
        [~, testInd] = min(abs(mFound - 15)); % Choose test peak index.
        testBack = trackData(:, testInd);
    catch
        disp('No test data found.');
    end
    testSub = testBack - (ArFilt .* (testBack(1) / ArFilt(1)));
    plot(times, testBack);
    hold on
    plot(times, testSub);
    box off;
    xlabel('Time [hours]');
    ylabel('Counts [a.u.]');
    legend('Raw', 'Baseline corr.');
    legend box off
catch
    disp('Normalisation error.');
end

%% Maybe come back to this.
% for iProcess = 1 : numel(mFound) % Iterate through peaks.
%     mLabel = find(ismember(mFound(iProcess), mRefs)); % Check reference.
%     if numel(mLabel) == 1
%         scaleData(:, iProcess) = scaleData(:, iProcess)
%     else
%     end
% end
% mTot = sum(trackData, 2); % Sum all masses in each spectrum.
