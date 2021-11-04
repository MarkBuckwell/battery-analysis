clear
clc
close all
fclose all;
addpath(cd);
set(groot, 'defaultFigureWindowStyle', 'Docked','defaultFigureColor', 'White',...
    'defaultAxesFontSize', 12, 'defaultAxesFontName', 'Arial',...
    'defaultLineLineWidth', 1.3, 'defaultAxesTickDir', 'Out',...
    'defaultAxesTickDirMode', 'Manual', 'defaultAxesLineWidth', 1.5,...
    'defaultFigureColormap', gray(64));

%% Mark Buckwell, Electrochemical Engineering Lab
% Chemical Engeineering Department, University College London
% March 2021, v1.0

%% Analysis of casing thickness for batches of reconstructed Nikon scans.
% Reads in 8-bit uint volumes, crops the top and bottom of the cell off,
% finds the peaks describing the casing for central X and Y cross-sections,
% gets the peak widths. So altogether this will give the distribution of
% casing thickness for each cell and for all cells input in the batch.
% Outlier data is removed as very high/low values for the casing width
% likely result from how Matlab defines the peak width, rather than real
% significant variations in the casing thickness.

%% Define parameters for processing.
% This will ideally be automated, either from file metadata or other source.
dataType = '*uint8'; % Data type to import as.
XPx = 727; % File pixel width in X.
YPx = 727; % File pixel width in Y.
ZPx = 2027; % File pixel width in Z.
topCrop = 400; % Slices to crop from top of data.
bottomCrop = 400; % Slices to crop from bottom of data.
voxelSize = 34.999; % Voxel dimension in um.
% Standard voxel size I use for P42As is 37.004 um.
% Standard voxel size I use for VTC6s is 34.999 um.
caseVal = 0.9; % Threshold normalised intensity to define peak as casing.
nSkip = 20; % Sampling interval in Z.
pStart = 240; % Row/column index to start peak finding.
pStop = 500; % Final row/column index for peak finding.
lPer = 7; % Lower percentile threshold to ignore outliers.
uPer = 70; % Upper percentile threshold to ignore outliers.
% 7/70 good for P42As.
% 7/70 good for VTC6s.

%% Choose file(s).
[fNames, fLoc] = uigetfile('*.*', 'MultiSelect', 'On'); % Read in files.
cd(fLoc); % Set directory to file location.
if ischar(fNames) == 1 % Put a single file in a cell.
    fNames = {fNames};
end

%% Process data.
fullMetrics = repmat({''}, numel(fNames), 5); % For all peak widths.
warning off
figure;
hold on
xlabel('Z slice');
ylabel('Casing thickness [um]');
ylim([3, 15] .* voxelSize);
fColours = hsv(numel(fNames));
for iFiles = 1 : numel(fNames) % Iterate through files.
    fullMetrics{iFiles, 1} = fNames{iFiles}; % Store file name.
    fID = fopen([fLoc, fNames{iFiles}]); % Open file.
    voxData = fread(fID, dataType); % Read in data.
    try % Reshape data into a volume.
        voxData = reshape(voxData, [XPx, YPx, ZPx]);
    catch % Adjust if input shape is wrong.
        voxData = reshape(voxData, [728, 728, 2026]);
        disp('Wrong shape, scanned for correct shape.');
    end
    % Crop top and bottom of cell and convert to single precision.
    voxData = single(voxData(:, :, topCrop : ZPx - bottomCrop));
	caseMetrics = zeros(ceil(size(voxData, 3) / nSkip), 4); % Array for final case widths.
	cleanMetrics = repmat({''}, 4, 1); % Ignoring outlier widths.
    for iZ = 1 : nSkip : size(voxData, 3) % Iterate through slices.
        rowMetrics = zeros(pStop - pStart + 1, 4);
        colMetrics = zeros(pStop - pStart + 1, 4);
        for iRows = pStart : pStop % Iterate through rows.
            [~, cLocs, cWidths] = findpeaks(...
                normalize(voxData(iRows, :, iZ), 'Range'),...
                'MinPeakHeight', caseVal, 'WidthReference', 'HalfHeight');
            if ~isempty(cWidths) == 1 % Check for presence of casing.
                [mVal, mInd] = min(cLocs); % Find first peak.
                rowMetrics(iRows - pStart + 1, 1) = cWidths(mInd); % Record width.
                rowMetrics(iRows - pStart + 1, 2) = mVal; % Record position.
                [mVal, mInd] = max(cLocs); % Find last peak.
                rowMetrics(iRows - pStart + 1, 3) = cWidths(mInd); % Record width.
                rowMetrics(iRows - pStart + 1, 4) = mVal; % Record position.
            end
        end
        % Get separation between most distant peaks.
        widthSeps = rowMetrics(:, 4) - rowMetrics(:, 2);
        [~, lInd] = max(widthSeps); % Find largest separation increasing.
        [~, uInd] = max(flipud(widthSeps)); % Find largest decreasing.
        caseMetrics(ceil(iZ / nSkip), 1 : 2) = rowMetrics(... % Store widths.
            ceil(numel(widthSeps) - uInd - (lInd / 2)), [1, 3]);
        for iCols = pStart : pStop % Iterate through columns.
            [~, cLocs, cWidths] = findpeaks(...
                normalize(voxData(:, iCols, iZ), 'Range'),...
                'MinPeakHeight', caseVal, 'WidthReference', 'HalfHeight');
            if ~isempty(cWidths) == 1 % Check for presence of casing.
                [mVal, mInd] = min(cLocs); % Find first peak.
                colMetrics(iCols - pStart + 1, 1) = cWidths(mInd); % Record width.
                colMetrics(iCols - pStart + 1, 2) = mVal; % Record position.
                [mVal, mInd] = max(cLocs); % Find last peak.
                colMetrics(iCols - pStart + 1, 3) = cWidths(mInd); % Record width.
                colMetrics(iCols - pStart + 1, 4) = mVal; % Record position.
            end
        end
        % Get separation between most distant peaks; closest to proper diameter.
        widthSeps = colMetrics(:, 4) - colMetrics(:, 2);
        [~, lInd] = max(widthSeps); % Find largest separation increasing.
        [~, uInd] = max(flipud(widthSeps)); % Find largest decreasing.
        % Take the midpoint between these.
        caseMetrics(ceil(iZ / nSkip), 3 : 4) = colMetrics(...  % Store widths.
            ceil(numel(widthSeps) - uInd - (lInd / 2)), [1, 3]);
        % Keep user updated on progress; file number, Z slice, widths.
        disp(strcat(string(iFiles), {', '}, string(iZ),...
            {', '}, string(caseMetrics(ceil(iZ / nSkip), 1)),...
            {', '}, string(caseMetrics(ceil(iZ / nSkip), 3))));
        % Plot the peak widths as they're found, to check on progress.
        plot(iZ, caseMetrics(ceil(iZ / nSkip), 1) .* voxelSize,...
            'o', 'Color', fColours(iFiles, :));
        plot(iZ, caseMetrics(ceil(iZ / nSkip), 3) .* voxelSize,...
            'o', 'Color', fColours(iFiles, :));
        pause(0.001);
    end
    fullMetrics{iFiles, 2} = caseMetrics; % Store raw data for file.
	% Replace outliers with NaN and store as scaled data.
    fullMetrics{iFiles, 3} = filloutliers(caseMetrics, nan,...
        'Percentiles', [lPer, uPer]) .* voxelSize;
    % Calculate mean and standard deviation of thickness for each file.
    fullMetrics{iFiles, 4} = mean(fullMetrics{iFiles, 3}, 'All', 'OmitNaN');
    fullMetrics{iFiles, 5} = std(fullMetrics{iFiles, 3}, 1, 'All', 'OmitNaN');
end

%% Plot all data to check.
figure;
subplot(2, 1, 1);
hold on
xlabel('Z slice');
ylabel('Casing thickenss [px]');
title('Raw data');
subplot(2, 1, 2);
hold on
xlabel('Z slice');
ylabel('Casing thickenss [um]');
title('Cleaned data');
for i = 1 : numel(fNames)
    subplot(2, 1, 1);
    plot(fullMetrics{i, 2}, 'o');
    subplot(2, 1, 2);
    plot(fullMetrics{i, 3}, 'o');
end

%% Save data as a table.
