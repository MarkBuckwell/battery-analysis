clear
clc
close all
addpath(cd);
set ( 0 , 'DefaultFigureWindowStyle' , 'Docked' ) ;
%% This function will unroll the jellyroll structure of cylindrical fuel cells.
% The user may load in as many cells of interest at once as desired. Each
% cell file should be a 2D image, in any format. The function will open the
% files sequentially and then find the centre of the cell. This is done by
% finding the widest point of each cell, i.e. the maximum difference
% between the peaks corresponding to the case. Note that the cell casing
% can have some width, as the function also finds the midpoint of the
% casing and accounts for this in defining the edges of the cell of the image.
% The function then takes the intersect between the widest points
% vertically and horizontally and defines this as the centre (this may be a
% simplification, something to develop later, or at a GUI for the user to
% itneract with the definition. The function then crops the image to one
% one pixel fewer than the determined width, and adds in order to have
% a central pixel to use as reference point. It then iterates circularly around
% the matrix holding the image data, finding the closest pixel to the
% required radial position, and plots the pixel value as a function of the
% angle and radius. The radial increment may be defined, but is as standard
% set to 1 (i.e. a single pixel).

%% Load in image files.
[FileGroup, DataPath] = uigetfile('*.*', 'DialogTitle',...
    'Select files:', 'MultiSelect', 'on'); % Gets file names and location.
cd(DataPath);
if double(ischar(FileGroup)) == 1
    NFC = 1;
else
    NFC = size(FileGroup, 2);
end
DataMatrix = repmat({''}, 1, NFC);
AllDiffs = repmat({''}, 2, NFC);
% Number of bright features to use in finding cell size, e.g. casing, mandrel.
% Too few features means the function might not find the edges. Only 2 are
% needed, but if the mandrel or ISC are white in the image, then these may
% be located instead of the casing. If more than 2 features are found, then
% the middle ones will be ignored in processing, assuming just the outer
% ones are relevant.
NFeatures = 8;
% Required feature height/intensity to classify as a data peak.
MinFeatureHeight = 0.99;
 for iNFC = 1 : NFC
	if NFC > 1 % For multiple files.
        FileChoice = strcat(DataPath, FileGroup{iNFC});
        FigTitle = FileGroup{iNFC}(1:end-4);
    else % For a single file.
        FileChoice = strcat(DataPath, FileGroup);
        FigTitle = FileGroup{iNFC}(1:end-4);
	end
	RollImage = imread(FileChoice);
	RollData = im2double(RollImage(:, :, 1));

    %% Find centre of image and edges of case.
    xDiffs = zeros(size(RollData, 1), 1);
    xWidths = zeros(size(RollData, 1), 2);
    yDiffs = zeros(size(RollData, 2), 1);
    yWidths = zeros(size(RollData, 2), 2);

    warning('Off', 'All'); % Ignore cases without peaks, i.e. outside roll.
    for i = 1 : size(RollData, 1) % Iterate through y, looking in x direction.
        % Find feature peaks for each line across the image.
        [~, Locs1] = findpeaks(RollData(i, :),...
            'NPeaks', NFeatures, 'MinPeakHeight', MinFeatureHeight);
        % Look in the other direction in order to find the feature width.
        [~, Locs2] = findpeaks(fliplr(RollData(i, :)),...
            'NPeaks', NFeatures, 'MinPeakHeight', MinFeatureHeight);
        % Flip the second feak find back round.
        Locs2= size(RollData, 2) + 1 - fliplr(Locs2);
        % In case fewer than 2 peaks are found.
        if numel(Locs1) && numel(Locs2) < 2
            xDiffs(i) = 0;
        elseif numel(Locs1) < 2 || numel(Locs2) < 2
            if numel(Locs2) > numel(Locs1)
                xDiffs(i) = diff(Locs2);
            elseif numel(Locs1) > numel(Locs2)
                xDiffs(i) = diff(Locs1);
            end
        else
            % Make sure only the outer features (casing) are used.
            Locs1 = [Locs1(1), Locs1(end)];
            Locs2 = [Locs2(1), Locs2(end)];
            xWidths(i, :) = Locs2 - Locs1;
            xDiffs(i) = diff(Locs2 - (xWidths(i, :) / 2));
        end
    end
    % Find centre point.
    [~, xCentre1] = max(xDiffs);
    [~, xCentre2] = max(flipud(xDiffs));
    xCentre = round(xCentre1 + (size(RollData, 2)...
        + 1 - xCentre2 - xCentre1) / 2);
    % Find edges of cell.
    xEdges = find(xDiffs);
    xEdges = [xEdges(1) - xWidths(xCentre, 1),...
        xEdges(end) - xWidths(xCentre, 2)];

    for i = 1 : size(RollData, 2) % Iterate through x, looking in y direction.
        [~, Locs1] = findpeaks(RollData(:, i),...
            'NPeaks', NFeatures, 'MinPeakHeight', MinFeatureHeight);
        [~, Locs2] = findpeaks(flipud(RollData(:, i)),...
            'NPeaks', NFeatures, 'MinPeakHeight', MinFeatureHeight);
        Locs2= size(RollData, 1) + 1 - flipud(Locs2);
        if numel(Locs1) && numel(Locs2) < 2
            yDiffs(i) = 0;
        elseif numel(Locs1) < 2 || numel(Locs2) < 2
            if  numel(Locs2) > numel(Locs1)
                yDiffs(i) = diff(Locs2);
            elseif numel(Locs1) > numel(Locs2)
                yDiffs(i) = diff(Locs1);
            end
        else
            Locs1 = [Locs1(1), Locs1(end)];
            Locs2 = [Locs2(1), Locs2(end)];
            yWidths(i, :) = Locs2 - Locs1;
            yDiffs(i) = diff(Locs2 - (yWidths(i, :) / 2));
        end
    end
    [~, yCentre1] = max(yDiffs);
    [~, yCentre2] = max(flipud(yDiffs));
    yCentre = round(yCentre1 + (size(RollData, 1)...
        + 1 - yCentre2 - yCentre1) / 2);
    yEdges = find(yDiffs);
    yEdges = [yEdges(1) - yWidths(yCentre, 1),...
        yEdges(end) - yWidths(yCentre, 2)];
    warning('On', 'All'); % Switch warnings back on.

    %% Crop image using narrower of x and y widths found.
    if diff(xEdges) > diff(yEdges)
        ImageSize = diff(yEdges);
    else
        ImageSize = diff(xEdges) - 1;
    end
    if mod(ImageSize, 2) > 0
        ImageSize = ImageSize - 1;
    end
    RollData(1:yCentre - (ImageSize / 2), :) = [];
    RollData(ImageSize:end, :) = [];
    RollData(:, 1:xCentre - (ImageSize / 2)) = [];
    RollData(:, ImageSize:end) = [];

    %% Unroll structure.
    rIncrement = 1; % Radius increment to use
    UnRollData = zeros(ImageSize / (2 * rIncrement), 360);
    r = 1;
    while r <= (ImageSize / 2) - 1 % Iterate increasing radius.
        for j = 1 : 360 % Iterate azimuth.
            x = (ImageSize / 2) + round(r * sind(j));
            y = (ImageSize / 2) + round(r * cosd(j));
            UnRollData(r, j) = RollData(y, x);
        end
        r = r + rIncrement;
    end
    figure(); % Plot!
    imagesc(UnRollData);
    xlabel('Angle [/theta]');
    ylabel('Radius [pixels]');
    set(gca, 'FontSize', 14, 'FontName', 'Helvetica',...
        'LineWidth', 1.2, 'Tickdir', 'Out', 'Box', 'Off');
    set(gcf, 'Color', 'w');
    title(FigTitle);
    DataMatrix{iNFC} = UnRollData; % Store the data for later.
    AllDiffs{1, iNFC} = xDiffs; % USeful for checking process.
    AllDiffs{2, iNFC} = yDiffs;
 end
