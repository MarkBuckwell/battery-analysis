function [CropCorrMatrix] = GaborCorrelationDetrender(iRot, iFit, iTrend)
%% Function to import Anand's csv files containins Gabor magnitude correlation data.
% The .csv files are pre-determined to a given format, which this script
% adheres to. Single or multiple files may be processed at once. For each
% file, the user selects a cropping region, implements smoothing if
% necessary, and ensures that the correct initial peak or valley is picked.
% They then confirm that the correct trend has been highlighted. For each
% file, a video may be produced produces and figures are saved for the
% correlation mapping, cropping, smoothing, data loss and failure
% propagation rate.

clc
%% Set default plotting and figure save options.
set(groot, 'defaultFigureWindowStyle', 'Docked','defaultFigureColor', 'White',...
    'defaultAxesFontSize', 14, 'defaultAxesFontName', 'Arial',...
    'defaultLineLineWidth', 2, 'defaultAxesTickDir', 'Out',...
    'defaultAxesTickDirMode', 'Manual', 'defaultAxesLineWidth', 2,...
    'defaultFigureColormap', jet(64));
SaveFormat = 'tiff'; % Set global figure save format.
addpath(cd);

%% Import and crop data from files.
% GaborImportCropFunction(SaveFormat), SaveFormat defines save format for
% figures. e.g. 'tiff', 'jpeg', 'eps'. This is defined globally above.
[CorrInfoMatrix, nFiles, FileGroup] = GaborCSVReader();
[CropCorrMatrix] = GaborCorrCropper(CorrInfoMatrix, nFiles,...
    FileGroup, SaveFormat);

%% Perform smoothing, fitting and trendfinding.
% iRot determines whether to find X in time or time in X:
% 1 - X in time;, 2 - time in X.
% iFit determines processing type:
% 1 - process data line-by-line; 2 - use image processing.
% iTrend determines means of finding trends.
% 1 - max/min values; 2 - max/min from multiple trends; 3 - findpeaks.
for iFitting = 1 : nFiles
    if iRot == 2 % Rotate data to look for time in X.
        CropCorrMatrix([1, 2]) = CropCorrMatrix([2, 1]);
        CropCorrMatrix{3} = CropCorrMatrix{3}';
        XLabel = 'X Cross-section (mm)';
        YLabel = 'Time (s)';
    else
        XLabel = 'Time (s)';
        YLabel = 'X Cross-section (mm)';
    end
    [CropCorrMatrix{iFitting, 4}, nLines] = GaborSmoother(...
        CropCorrMatrix(iFitting, 1 : 3), FileGroup, nFiles, iFitting,...
        SaveFormat, YLabel); % Perform data smoothing.
    % User chooses to search for a valley or peak trend.
    iPeak = 0;
    TrendChoice = questdlg('Search for peak or valley trend?',...
        'Trend selection', 'Peak trend/rising edge',...
        'Valley trend/falling edge', 'Peak trend/rising edge');
    switch TrendChoice
        case 'Peak trend/rising edge'
            iPeak = 1; % Do nothing to data.
        case 'Valley trend/falling edge'
            iPeak = -1; % Invert data.
    end
    
%% Choose whether or not to capture video.
    VideoChoice = questdlg('Capture trend find videos?', 'Video selection',...
        'Yes', 'No', 'No');
    switch VideoChoice
        case 'Yes'
            if nFiles == 1 % Prepare a video.
                PeakVideo = VideoWriter('Trend_Video', 'MPEG-4');
            else
                PeakVideo = VideoWriter(strcat(FileGroup{iFitting},...
                    '_Trend_Video'), 'MPEG-4');
            end
            open(PeakVideo);
            iVideo = 1;
        case 'No' % Pass something empty.
            PeakVideo = {''};
            iVideo = 0;
    end
    
%% Find data trend(s);
    if iFit == 1 % Fit line-by-line.
        if iTrend == 1 % Max/min X value at each timepoint.
            CropCorrMatrix{iFitting, 5} = GaborMaxMin(...
                CropCorrMatrix(iFitting, :), FileGroup, nFiles,...
                nLines, iPeak, PeakVideo, iVideo, iFitting, SaveFormat,...
                XLabel, YLabel);
        elseif iTrend == 2 % Max/min from multiple trends at each timepoint.
            CropCorrMatrix{iFitting, 5} = GaborMultiMaxMin(...
                CropCorrMatrix(iFitting, :), FileGroup, nFiles,...
                nLines, iPeak, PeakVideo, iVideo, iFitting, SaveFormat,...
                XLabel, YLabel);
        elseif iTrend == 3 % Peak finding per timepoint.
            CropCorrMatrix{iFitting, 5} = GaborPeakFind(...
                CropCorrMatrix(iFitting, :), FileGroup, nFiles,...
                nLines, iPeak, PeakVideo, iVideo, iFitting, SaveFormat,...
                XLabel, YLabel);
        end
    elseif iFit == 2 % Fit using image processing.
        return % Do nothing for now
    end
    if iVideo == 1 % Save video.
        close(PeakVideo)
    end
    
%% Plot and fit data trend(s).
    CropCorrMatrix{iFitting, 6} = GaborTrendPlotter(...
        CropCorrMatrix(iFitting, [1, 5]), FileGroup, nFiles, iFitting,...
        SaveFormat, XLabel, YLabel, iRot);
    
%% Save processed data.
    if nFiles == 1
        save('Processed_Data.mat', 'CropCorrMatrix');
    else
        save(strcat(FileGroup{iFitting}, '_Processed_Data.mat'),...
            'CropCorrMatrix');
    end
end
end

%% GaborCSVReader imports files.
function [CorrInfoMatrix, nFiles, FileGroup] = GaborCSVReader()
[FileGroup, DataPath] = uigetfile('*.csv', 'DialogTitle',...
    'Select files:', 'MultiSelect', 'on'); % Gets file names and location.
cd(DataPath);
SingleFile = double(ischar(FileGroup));
if  SingleFile > 0
    nFiles = 1;
else
    nFiles = length(FileGroup);  % Number of files to import.
end
CorrInfoMatrix = repmat({''}, nFiles, 2);
% Open files and get data.
if nFiles == 1
    FileID = fopen(FileGroup);
    HeaderLine = textscan(FileID, '%s', 'Delimiter', '\n');
    HeaderLine = textscan(HeaderLine{1}{1}, '%s', 'Delimiter', ',');
    CorrInfoMatrix{1} = HeaderLine{1}';
    CorrInfoMatrix{2} = readmatrix(FileGroup, 'NumHeaderLines', 1);
else
    for iRead = 1 : nFiles
        FileID = fopen(FileGroup{iRead});
        HeaderLine = textscan(FileID, '%s', 'Delimiter', '\n');
        HeaderLine = textscan(HeaderLine{1}{1}, '%s', 'Delimiter', ',');
        CorrInfoMatrix{iRead, 1} = HeaderLine{1}';
        CorrInfoMatrix{iRead, 2} = readmatrix(FileGroup{iRead}, 'NumHeaderLines', 1);
    end
end
end

%% GaborCorrCropper displays data, allows user to crop, then parses cropped data.
function [CropCorrMatrix] = GaborCorrCropper(CorrInfoMatrix, nFiles,...
    FileGroup, SaveFormat)
% Find gap markers in csv.
CropCorrMatrix = repmat({''}, nFiles, 6);
% 1st column - cropped time values; 2nd column - cropped X values; 3rd
% column - cropped correlation data, rotated by 90 degrees, and mirrored,
% if applicable; 4th column - smoothed data, 5the column - trend data.
for iCrop = 1 : nFiles
    GapMarkers = find(isnan(CorrInfoMatrix{iCrop, 2}(1, :)));
    CrossDataStart = GapMarkers(2) + 1;
    CrossDataEnd = GapMarkers(3) - 1;
    CorrData = rot90(CorrInfoMatrix{iCrop, 2}(:,...
        CrossDataStart : CrossDataEnd));
    CorrData = CorrData ./ max(max(CorrData)); % Normalise data.
    % Find measurement parameters and set minima to 0.
    tColumn = find(contains(CorrInfoMatrix{iCrop, 1}, 'Time (s)'));
    if isempty(tColumn) == 1 % Arbitrary time indices.
        tVals = 0 : size(CorrData, 2) - 1;
    else % Real time values.
        tVals = CorrInfoMatrix{iCrop, 2}(:, tColumn);
%         tVals = tVals - tVals(1);
    end
    XVals = textscan(cell2mat(CorrInfoMatrix{iCrop, 1}(CrossDataStart :...
        CrossDataEnd)), '%f', 'Delimiter', '""'); % X position values in mm.
    XVals = (rmmissing(XVals{1}) - min(XVals{1})) ./ 1000 ;
    % Show map and allow user to define crop region.
    % Get all data and show initial map.
    figure(1);
    clf
    xlabel('Time (s)');
    ylabel('X Cross-section (mm)');
    hold on
    imagesc([tVals(1), tVals(end)],...
        [XVals(1), XVals(end)], CorrData);
    xlim([tVals(1), tVals(end)]);
    ylim([XVals(1), XVals(end)]);
    cBar = colorbar;
    ylabel(cBar, 'Normalised GaborCorr Magnitude',...
        'FontSize', get(0, 'defaultAxesFontSize') - 2);
    hold off
    iMirror = 0;
    while iMirror < 1
        AskMirror = questdlg('Reverse X Cross-section?',...
            'Data reverser', 'Yes', 'No', 'Exit', 'No');
        switch AskMirror % Allow user to reverse X data.
            case 'Yes'
                CorrData = flipud(CorrData);
                clf
                xlabel('Time (s)');
                ylabel('X Cross-section (mm)');
                MyPlotOptions(get(0, 'defaultAxesFontSize'));                
                hold on
                imagesc([tVals(1), tVals(end)],...
                    [XVals(1), XVals(end)], CorrData);
                xlim([tVals(1), tVals(end)]);
                ylim([XVals(1), XVals(end)]);                    
                cBar = colorbar;
                ylabel(cBar, 'Normalised GaborCorr Magnitude',...
                    'FontSize', get(0, 'defaultAxesFontSize') - 2);
                hold off
            case 'No'
                iMirror = 1;
            case 'Exit'
                return
        end
    end
    if nFiles == 1 % Save figure.
        saveas(gcf, strcat('Full_Map.', SaveFormat), SaveFormat);
    else
        saveas(gcf, strcat(FileGroup{iCrop},...
            '_Full_Map.', SaveFormat), SaveFormat);
    end
    % Allow user to define crop region.
    iDefine = 0;
    % Pre-define cropping boundaries in case no cropping done.
    XMin = 1;
    XMax = numel(XVals);
    tMin = 1;
    tMax = numel(tVals);
    CropChoices = {num2str(XVals(1)) num2str(XVals(end))...
        num2str(tVals(1)) num2str(tVals(end))};
    while iDefine < 1
        AskCrop = questdlg('Accept or change cropping?',...
            'Crop selection', 'Crop data', 'Show original data',...
            'Accept', 'Show original data');
        switch AskCrop
            case 'Crop data' % Allow user to define cropping region.
                UserPrompt = {'X start position', 'X end position',...
                    'Start time', 'End time'};
                DefaultAnswers = {CropChoices{1} CropChoices{2}...
                    CropChoices{3} CropChoices{4}};
                CropChoices = inputdlg(UserPrompt,...
                    'Choose cropping paremeters', 1, DefaultAnswers);
                [~, XMin] = min(abs(XVals -...
                    str2double(CropChoices{1})));
                [~, XMax] = min(abs(XVals -...
                    str2double(CropChoices{2})));
                [~, tMin] = min(abs(tVals -...
                    str2double(CropChoices{3})));
                [~, tMax] = min(abs(tVals -...
                    str2double(CropChoices{4})));
                figure(2); % Show cropped region.
                clf
                xlabel('Time (s)');
                ylabel('X Cross-section (mm)');
                MyPlotOptions(get(0, 'defaultAxesFontSize'));
                hold on
                imagesc([tVals(tMin), tVals(tMax)],...
                    [XVals(XMin), XVals(XMax)],...
                    CorrData(XMin : XMax, tMin : tMax));
                xlim([tVals(tMin), tVals(tMax)]);
                ylim([XVals(XMin), XVals(XMax)]);
                cBar = colorbar;
                ylabel(cBar, 'Normalised GaborCorr Magnitude',...
                    'FontSize', get(0, 'defaultAxesFontSize') - 2);
            case 'Show original data' % Return to display original data.
                figure(1);
            case 'Accept'
                CropCorrMatrix{iCrop, 1} = tVals(tMin : tMax);
                CropCorrMatrix{iCrop, 2} = XVals(XMin : XMax);
                CropCorrMatrix{iCrop, 3} = CorrData(XMin : XMax,...
                    tMin : tMax);
                iDefine = 1;
                % Save figure(s).
                if nFiles == 1
                    saveas(gcf, strcat('Cropped_Map.', SaveFormat),...
                        SaveFormat);
                else
                    saveas(gcf, strcat(FileGroup{iCrop},...
                        '_Cropped_Map.', SaveFormat), SaveFormat);
                end
        end
    end
end
end

%% GaborSmoother allows user to define data smoothing parameters.
function [SmoothedData, nLines] = GaborSmoother(CropCorrData, FileGroup,...
    nFiles, iFitting, SaveFormat, YLabel)
nLines = numel(CropCorrData{1});
figure(3);
clf
LineOut = CropCorrData{3}(:, ceil(nLines / 2));
LinePlot = plot(CropCorrData{2}, LineOut);
set(LinePlot, 'Color', 'Black', 'LineStyle', '-');
MyPlotOptions(get(0, 'defaultAxesFontSize'));
xlabel(YLabel);
ylabel('Normalised GaborCorr Magnitude');
hold on
iSmooth = 0 ;
while iSmooth < 1
    LineChoice = questdlg('Smooth data?', 'Smoothing choice',...
        'Yes, use this line', 'Look at another line', 'No',...
        'Look at another line');
    switch LineChoice
        case 'Yes, use this line'
            SmoothChoices = {'3', '3'}; % Initial smoothing parameters.
            iFilter = 0;
            while iFilter < 1
                SmoothChoices = inputdlg(...
                    {'1) lowpass, 2) median, 3) sgolay',...
                    'Filter width'}, 'Select filtering parameters',...
                    1, SmoothChoices);
                try
                    if str2double(SmoothChoices{1}) == 1
                        SmoothLine = lowpass(LineOut,...
                            str2double(SmoothChoices{2}));
                    elseif str2double(SmoothChoices{1}) == 2
                        SmoothLine = medfilt1(LineOut,...
                            str2double(SmoothChoices{2}));
                    elseif str2double(SmoothChoices{1}) == 3
                        SmoothLine = smooth(LineOut, 'sgolay',...
                            str2double(SmoothChoices{2}));
                    end
                    SmoothPlot = plot(CropCorrData{2}, SmoothLine);
                    set(SmoothPlot, 'Color', 'Red', 'LineStyle', '-');
                catch
                    disp('Invalid filter width, please try again.');
                end
                AskProceed = questdlg('Proceed with smoothing data?',...
                    'Smooth check', 'Yes', 'Change parameters',...
                    'No, exit', 'Change parameters');
                switch AskProceed
                    case 'Yes'
                        if nFiles == 1 % Save figure.
                            saveas(gcf, strcat('Smooth_Line.',...
                                SaveFormat), SaveFormat);
                        else
                            saveas(gcf, strcat(FileGroup{iFitting},...
                                '_Smooth_Line.', SaveFormat), SaveFormat);
                        end
                        try % In case fitting accepting with incorrect parameters.
%                         delete(LinePlot);
%                         delete(SmoothPlot);
                        catch % Dont't do anything.
                        end
                        SmoothedData = zeros(numel(CropCorrData{2}),...
                            nLines);
%                         ylabel('Data Loss (data - smoothed)');
                        for iLines = 1 : nLines
                            LineOut = CropCorrData{3}(:, iLines);
                            if str2double(SmoothChoices{1}) == 1
                                SmoothLine = lowpass(LineOut,...
                                    str2double(SmoothChoices{2}));
                            elseif str2double(SmoothChoices{1}) == 2
                                SmoothLine = medfilt1(LineOut,...
                                    str2double(SmoothChoices{2}));
                            elseif str2double(SmoothChoices{1}) == 3
                                SmoothLine = smooth(LineOut, 'sgolay',...
                                    str2double(SmoothChoices{2}));
                            end
                            SmoothedData(:, iLines) = SmoothLine;
%                             plot(CropCorrData{2},... % Plot residual.
%                                 CropCorrData{3}(:, iLines) - SmoothLine);
%                             pause(0.01);
                        end
%                         if nFiles == 1 % Save figure.
%                             saveas(gcf, strcat('Smooth_All.',...
%                                 SaveFormat), SaveFormat);
%                             saveas(gcf, 'Smooth_All.fig', 'fig');
%                         else
%                             saveas(gcf, strcat(FileGroup{iFitting},...
%                                 '_Smooth_All.', SaveFormat), SaveFormat);
%                             saveas(gcf, strcat(FileGroup{iFitting},...
%                                 '_Smooth_All.fig'), 'fig');
%                         end
                        iSmooth = 1;
                        iFilter = 1;
                    case 'Change parameters'
                        try % In case fitting accepted with incorrect parameters.
                            delete(SmoothPlot);
                        catch % Don't do anything if SmoothPlot doesn't exist.
                        end
                    case 'No, exit'
                        SmoothedData = CropCorrData{3};
                        return
                end
            end
        case 'Look at another line'
            LineOut = CropCorrData{3}(:, randi(nLines));
            delete(LinePlot);
            LinePlot = plot(CropCorrData{2}, LineOut);
        case 'No'
            SmoothedData = CropCorrData{3};
            iSmooth = 1;
    end
end
end

%% GaborMaxMin identifies trends using max/min values.
function [PeakInfo] = GaborMaxMin(CropCorrData, FileGroup, nFiles,...
    nLines, iPeak, PeakVideo, iVideo, iFitting, SaveFormat, XLabel, YLabel)
PeakInfo = zeros(nLines, 1);
figure(4);
hold on
imagesc([CropCorrData{1}(1), CropCorrData{1}(end)],...
    [CropCorrData{2}(1), CropCorrData{2}(end)], CropCorrData{3});
xlim([CropCorrData{1}(1), CropCorrData{1}(end)]);
ylim([CropCorrData{2}(1), CropCorrData{2}(end)]);
ylabel(YLabel);
xlabel(XLabel);
cBar = colorbar;
ylabel(cBar, 'Normalised Correlation Magnitude',...
    'FontSize', get(0, 'defaultAxesFontSize') - 2);
for iFind = 1 : nLines
    LineOut = CropCorrData{4}(:, iFind) .* iPeak;
    [~, MInd] = max(LineOut);
    PeakInfo(iFind) = CropCorrData{2}(MInd);
    TrendPoint = plot(CropCorrData{1}(iFind), PeakInfo(iFind), 'o');
    if iFind == 1
        set(TrendPoint, 'Color', 'Black');
    else
        set(TrendPoint, 'Color', 'White');
    end
    pause(0.01);
    if iVideo == 1
        writeVideo(PeakVideo, getframe(gcf));
    end
end
if nFiles == 1 % Save figure.
    saveas(gcf, strcat('Trend_Fit_Map.', SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iFitting}, 'Trend_Fit_Map.', SaveFormat),...
        SaveFormat);
end
end

%% GaborMultiMaxMin identifies trends using max/min values where multiple
% trends are present. Currently, only a single trend is selected.
function [PeakInfo] = GaborMultiMaxMin(CropCorrData, FileGroup, nFiles,...
    nLines, iPeak, PeakVideo, iVideo, iFitting, SaveFormat, XLabel, YLabel)
PeakInfo = zeros(nLines, 1);
figure(4);
hold on
imagesc([CropCorrData{1}(1), CropCorrData{1}(end)],...
    [CropCorrData{2}(1), CropCorrData{2}(end)], CropCorrData{3});
xlim([CropCorrData{1}(1), CropCorrData{1}(end)]);
ylim([CropCorrData{2}(1), CropCorrData{2}(end)]);
ylabel(YLabel);
xlabel(XLabel);
cBar = colorbar;
ylabel(cBar, 'Normalised Correlation Magnitude',...
    'FontSize', get(0, 'defaultAxesFontSize') - 2);
MyPlotOptions(get(0, 'defaultAxesFontSize'));
PeakChoices = {'0.5' '2'}; % Peak find parameters. 1 - starting position of
% peak; 2 - range to search in.
iCheck = 0;
while iCheck < 1
    PeakChoices = inputdlg({'Peak starting position in X/time',...
        'Search range'}, 'Choose peak finding parameters', 1, PeakChoices);
	[~, PRange] = min(abs(CropCorrData{2}...
        - (str2double(PeakChoices{2}) + CropCorrData{2}(1))));
    % Range of matrix elements to check.
    PRange = floor(PRange / 2);
    PPos = str2double(PeakChoices{1});
    for iFind = 1 : nLines
        LineOut = CropCorrData{4}(:, iFind) .* iPeak;
        [~, PPos] = min(abs(CropCorrData{2} - PPos));
        LowerBound = PPos - PRange;
        if LowerBound < 1 % Ensure range doesn't go below 1st index.
            LowerBound = 1;
        end
        UpperBound = PPos + PRange; % Ensure range doesn't go beyond final index.
        if UpperBound > numel(CropCorrData{2})
            UpperBound = numel(CropCorrData{2});
        end            
        [~, MInd] = max(LineOut(LowerBound : UpperBound));
        PeakInfo(iFind) = CropCorrData{2}(MInd + LowerBound - 1);
        PPos = PeakInfo(iFind); % Retain previous value.
        TrendPoint = plot(CropCorrData{1}(iFind), PeakInfo(iFind), 'o');
        if iFind == 1
            set(TrendPoint, 'Color', 'Black');
        else
            set(TrendPoint, 'Color', 'White');
        end
        pause(0.01);
        if iVideo == 1
            writeVideo(PeakVideo, getframe(gcf));
        end
    end
    AcceptTrend = questdlg('Accept located trend?', 'Trend check', 'Yes',...
            'No, change parameters', 'Exit', 'No, change parameters');
    switch AcceptTrend
        case 'Yes'
            iCheck = 1;
        case 'No, change parameters'
            clf
            hold on
            imagesc([CropCorrData{1}(1), CropCorrData{1}(end)],...
                [CropCorrData{2}(1), CropCorrData{2}(end)], CropCorrData{3});
            xlim([CropCorrData{1}(1), CropCorrData{1}(end)]);
            ylim([CropCorrData{2}(1), CropCorrData{2}(end)]);
            ylabel(YLabel);
            xlabel(XLabel);
            cBar = colorbar;
            ylabel(cBar, 'Normalised Correlation Magnitude',...
                'FontSize', get(0, 'defaultAxesFontSize') - 2);
            MyPlotOptions(get(0, 'defaultAxesFontSize'));
    end
end
if nFiles == 1 % Save figure.
    saveas(gcf, strcat('Trend_Fit_Map.', SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iFitting}, 'Trend_Fit_Map.', SaveFormat),...
        SaveFormat);
end
end

%% GaborPeakFind identifies trends using findpeaks.
function [PeakInfo] = GaborPeakFind(CropCorrData, FileGroup, nFiles, nLines,...
    iPeak, PeakVideo, iVideo, iFitting, SaveFormat, XLabel, YLabel)
PeakInfo = zeros(nLines, 1);
figure(4);
hold on
imagesc([CropCorrData{1}(1), CropCorrData{1}(end)],...
    [CropCorrData{2}(1), CropCorrData{2}(end)], CropCorrData{3} .* iPeak);
xlim([CropCorrData{1}(1), CropCorrData{1}(end)]);
ylim([CropCorrData{2}(1), CropCorrData{2}(end)]);
ylabel(YLabel);
xlabel(XLabel);
cBar = colorbar;
ylabel(cBar, 'Normalised Correlation Magnitude',...
    'FontSize', get(0, 'defaultAxesFontSize') - 2);
PeakChoices = {'2' '1' '0.5'}; % Peak find parameters. 1 - number of peaks
% to find, 2 - chosen peak to start trend, 3 - minimum peak height.
iCheck = 0;
while iCheck < 1
    PeakChoices = inputdlg({'Number of peaks to find', 'Chosen peak',...
        'Minimum peak height'}, 'Choose peak finding parameters', 1,...
        PeakChoices);
    try
        TrendPoints = repmat({''}, nLines, 1);
        for iFind = 1 : nLines
            LineOut = CropCorrData{4}(:, iFind) .* iPeak;
            [~, PInd] = findpeaks(LineOut,'NPeaks', str2double(PeakChoices{1}),...
                    'MinPeakHeight', str2double(PeakChoices{3}));
            if PInd(1) ~= 0
                if  iFind == 1 % Use chosen peak.
                    Closest = PInd(str2double(PeakChoices{2}));
                    PeakInfo(iFind) = CropCorrData{2}(Closest);
                else % Find closest peak to previous.
                    [~, Closest] = min(abs(PInd - Closest));
                    PeakInfo(iFind) = CropCorrData{2}(Closest);
                end
                TrendPoints{iFind} = plot(CropCorrData{1}(iFind),...
                    PeakInfo(iFind), 'o');
                if iFind == 1
                    set(TrendPoints{iFind}, 'Color', 'Black');
                else
                    set(TrendPoints{iFind}, 'Color', 'White');
                end
                pause(0.001);
            end
        end
        TrendCheck = questdlg('Accept located trend?', 'Trend check', 'Yes',...
            'No, change parameters', 'Exit', 'No, change parameters');
        switch TrendCheck
            case 'Accept'
                iCheck = 1;
            case 'No, change parameters'
                delete(TrendPoints);
            case 'Exit'
                return
        end
    catch
        disp('Invalid search parameters, please try again.');       
    end
end
if iVideo == 1
    writeVideo(PeakVideo, getframe(gcf));
end
if nFiles == 1 % Save figure.
    saveas(gcf, strcat('Trend_Fit_Map.', SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iFitting}, 'Trend_Fit_Map.', SaveFormat),...
        SaveFormat);
end
end

%% GaborTrendPlotter plots final failure propagation trend(s).
function [PropagationData] = GaborTrendPlotter(CropCorrData, FileGroup,...
    nFiles, iFitting, SaveFormat, XLabel, YLabel, iRot)
%% Plot failure propagation data.
if iRot == 1 % For finding X in time.
    XVals = CropCorrData{2};
    tVals = CropCorrData{1};
else % Switch the data types if finding time in X.
    XVals = CropCorrData{1};
    tVals = CropCorrData{2};
end
XOffset = 0; % Define zero offset in all cases.
tOffset = 0;
AskZero = questdlg('Zero time/X?', 'For pinned fitting', 'Yes', 'No', 'Yes');
switch AskZero % Store offset and set intercept to 0 for fitting.
    case 'Yes'
        XOffset = XVals(1);
        XVals = XVals - XOffset;
        tOffset = tVals(1);
        tVals = tVals - tOffset;
end
iPlot = 0;
CropChoices = {num2str(tVals(1)) num2str(tVals(end))};
FitChoices = {'poly2', '-100, -100, -100', '100, 100, 100'};
while iPlot < 1
    figure(5);
    clf
    plot(tVals, XVals, 'o');
    ylabel('X Cross-section (mm)');
    xlabel('Time (s)');
    MyPlotOptions(get(0, 'defaultAxesFontSize'));
    if nFiles == 1 % Save figure.
        saveas(gcf, strcat('Full_Propagation_Plot.', SaveFormat), SaveFormat);
    else
        saveas(gcf, strcat(FileGroup{iFitting}, '_Full_Propagation_Plot.',...
            SaveFormat), SaveFormat);
    end

    %% Allow user to determine fitting region and options.
    iDefine = 0;
    tMin = 1;
    tMax = numel(tVals);
    while iDefine < 1
        AskCrop = questdlg('Accept or change cropping?',...
            'Crop selection', 'Crop data', 'Show original data',...
            'Accept', 'Show original data');
        switch AskCrop
            case 'Crop data' % Allow user to define cropping region.
                DefaultAnswers = {CropChoices{1} CropChoices{2}};
                CropChoices = inputdlg({'Start time/X', 'End time/X'},...
                    'Choose cropping paremeters', 1, DefaultAnswers);
                [~, tMin] = min(abs(tVals -...
                    str2double(CropChoices{1})));
                [~, tMax] = min(abs(tVals -...
                    str2double(CropChoices{2})));
                figure(6); % Show cropped region.
                clf
                plot(tVals(tMin : tMax), XVals(tMin : tMax), 'o');
                xlabel('Time (s)');
                ylabel('X Cross-section (mm)');
                MyPlotOptions(get(0, 'defaultAxesFontSize'));
            case 'Show original data' % Return to display original data.
                figure(5);
            case 'Accept'
                hold on
                clear PropagationData
                PropagationData(:, 1) = tVals(tMin : tMax);
                PropagationData(:, 2) = XVals(tMin : tMax);
                switch AskZero
                    case 'Yes'
                        XOffset = XOffset + XVals(tMin);
                        tOffset = tOffset + tVals(tMin);
                        PropagationData(:, 1) = PropagationData(:, 1)...
                            - tVals(tMin);
                        PropagationData(:, 2) = PropagationData(:, 2)...
                            - XVals(tMin);
                end
                iDefine = 1;
        end
    end

    %% Fit data.
    iFit = 0;
    while iFit < 1
        FitChoices = inputdlg({'Fit type', 'Lower Bound', 'Upper Bound'},...
            'Choose fitting paremeters, comma-separated', 1, FitChoices);
        % LAR seems most appropriate, as all data, points equally valid, no outliers.
        TFitOptions = fitoptions('Method', 'LinearLeastSquares',...
            'Robust', 'LAR',...
            'Lower', cell2mat(textscan(FitChoices{2}, '%f', 'Delimiter', ','))',...
            'Upper', cell2mat(textscan(FitChoices{3}, '%f', 'Delimiter', ','))');
        [TrendFit, GoF, FitOut] = fit(PropagationData(:, 1),...
            PropagationData(:, 2), FitChoices{1}, TFitOptions);
        figure(6);
        clf % Clear and replot with non-zeroed cropped data.
        plot(PropagationData(:, 1) + tOffset,...
            PropagationData(:, 2) + XOffset, 'o');
        hold on
        xlabel('Time (s)');
        ylabel('X Cross-section (mm)');
        MyPlotOptions(get(0, 'defaultAxesFontSize'));        
        FitPlot = plot(PropagationData(:, 1) + tOffset,...
            TrendFit(PropagationData(:, 1)) + XOffset);
        delete(legend);
        xlabel('Time (s)');
        ylabel('X Cross-section (mm)');
        AcceptFit = questdlg('Accept or change fitting parameters?',...
            'Fitting confirmation', 'Accept', 'Change parameters', 'Exit',...
            'Change parameters');
        switch AcceptFit
            case 'Accept'
                iFit = 1;
                iPlot = 1;
                PropagationData(:, 1) = PropagationData(:, 1) + tOffset;
                PropagationData(:, 2) = PropagationData(:, 2) + XOffset;
                PropagationData(:, 3) = FitOut.residuals;
                if nFiles == 1 % Save figure(s).
                    saveas(gcf, strcat('Cropped_Propagations_Plot.',...
                        SaveFormat), SaveFormat);
                    save('Fit_Out.mat', 'FitOut');
                    save('Trend_Fit.mat', 'TrendFit');
                    save('Goodness_of_Fit.mat', 'GoF');
                else
                    saveas(gcf, strcat(FileGroup{iCrop},...
                        '_Cropped_Propagation_Plot.', SaveFormat), SaveFormat);
                    save(strcat(FileGroup{iCrop}, '_Fit_Out.mat'), 'FitOut');
                    save(strcat(FileGroup{iCrop}, '_Trend_Fit.mat'), 'TrendFit');
                    save(strcat(FileGroup{iCrop}, '_Goodness_of_Fit.mat'), 'GoF');
                end
            case 'Change parameters'
                iFit = 1;
                delete(FitPlot);
            case 'Exit'
                return
        end
    end
end
% Overlay fit on previous plots.
figure(2);
plot(PropagationData(:, 1), TrendFit(PropagationData(:, 1)), 'White');
delete(legend);
xlabel('Time (s)');
ylabel('X Cross-section (mm)');
if nFiles == 1 % Save figure(s).
    saveas(gcf, strcat('Fit_on_Map.',...
        SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iCrop},...
        '_Fit_on_Map.', SaveFormat), SaveFormat);
end
figure(4);
if iRot == 2
    plot(TrendFit(tVals(tMin) : 0.001 : tVals(tMax)) + XOffset,...
        (tVals(tMin) : 0.001 : tVals(tMax)) + tOffset, 'Black');
else
    plot((tVals(tMin) : 0.001 : tVals(tMax)) + tOffset,...
        TrendFit(tVals(tMin) : 0.001 : tVals(tMax)) + XOffset, 'Black');    
end
delete(legend);
xlabel(XLabel);
ylabel(YLabel);
if nFiles == 1 % Save figure(s).
    saveas(gcf, strcat('Data_Fit_on_Map.',...
        SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iCrop},...
        '_Data_Fit_on_Map.', SaveFormat), SaveFormat);
end
figure(6); % Show differential of fit.
clf
yyaxis left
plot(PropagationData(:, 1), PropagationData(:, 2), 'o');
hold on
plot((tVals(tMin) : 0.001 : tVals(tMax)) + tOffset,...
    TrendFit(tVals(tMin) : 0.001 : tVals(tMax)) + XOffset, 'Black');
ylabel('X Cross-section (mm)');
yyaxis right
plot((tVals(tMin) : 0.001 : tVals(tMax)) + tOffset,...
    differentiate(TrendFit, tVals(tMin) : 0.001 : tVals(tMax)) + XOffset);
ylabel('Propagation rate (mms^{-1})');
xlabel('Time (s)');
delete(legend);
MyPlotOptions(get(0, 'defaultAxesFontSize'));
if nFiles == 1 % Save figure.
    saveas(gcf, strcat('Propagation_Plot.', SaveFormat), SaveFormat);
else
    saveas(gcf, strcat(FileGroup{iFitting}, '_Propagation_Plot.', SaveFormat),...
        SaveFormat);
end
end