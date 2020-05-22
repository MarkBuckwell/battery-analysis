clear
clc
close all
addpath(cd);
set (0, 'DefaultFigureWindowStyle', 'Docked');

%% This script takes a video of a layered cell undergoing failure, e.g. from
% high-speed x-ray imaging of thermal or mechanical abuse, and filters the
% images to better study the changing layers or remove, for example, the
% impinging nail. All images in the video are cropped based on the first
% frame, to remove any black areas surrounding the region of interest, plus
% an additional crop to remove any unwanted edges. Videos are produced of
% the filtering in progress.

%% Global plotting display parameters (currently unused).
% FigFontSize = 14;
% FigLineWidth = 1.2;
% FigFont = 'Helvetica';

%% Get video address and read in video.
[FileGroup, DataPath] = uigetfile('*.*', 'DialogTitle',...
    'Select files:', 'MultiSelect', 'off'); % Gets file names and location.
cd(DataPath);
VidRead = VideoReader(strcat(DataPath, FileGroup));
VidStep = 10; % Number of frames to step through per iteration.
%     LBPSurf = zeros(26, VidRead.NumFrames);
GFilter = 8; % Gabor filter size.
GStep = 10; % Angular step in Gabor filtering.
% Gabor filtered data.
GArray = zeros((180 / GStep) + 1, VidRead.NumFrames / VidStep);

%% Determine video cropping parameters.
CropFrame = read(VidRead, 1); % First frame, to determine cropping.
xTrim = 50; % Additional pixels to crop horixontally.
xCrop = CropFrame(size(CropFrame, 1) / 2, :, 1); % Line across middle of frame.
xCrop = find(xCrop); % Nonzero elements.
xCrop = [xCrop(1) + xTrim, xCrop(end) - xTrim]; % Indices to crop to.
yTrim = 50; % Additional pixels to crop vertically.
yCrop = CropFrame(:, size(CropFrame, 2) / 2, 1); % Line across middle of frame.
yCrop = find(yCrop); % Nonzero elements.
yCrop = [yCrop(1) + yTrim, yCrop(end) - yTrim]; % Indices to crop to.

%% Set up video outputs.
VidQuad = zeros(yCrop(end) - yCrop(1) + 1, xCrop(end) - xCrop(1) + 1, 1);
QuadVideo = VideoWriter('QuadVideo.avi');
QuadVideo.FrameRate = 2; % Frames per second.
open(QuadVideo);
GaborVideo = VideoWriter('GaborVideo.avi');
GaborVideo.FrameRate = 20; % Frames per second.
open(GaborVideo);

%% Process video.
for i = 0:VidStep:VidRead.NumFrames - 1 % Step through video frames.
    VidFrame = read(VidRead, i + 1); % Read frame.
    VidFrame = VidFrame(yCrop(1):yCrop(end), xCrop(1):xCrop(end), :);
    figure(1); % Graphic for raw video and Gabor filtering.
    subplot(2, 2, 1);
    imagesc(VidFrame); % Show video frame.
    title('Raw video');
    set(gca, 'XTick', [], 'YTick', []); % Remove axes and labels.
    
    %% Remove background from frame.
    VidRaw = double(VidFrame(:, :, 1));
    for j = 1:size(VidRaw, 2) % Step through columns of video frame.
        VidQuad(:, j) = detrend(VidRaw(:, j), 2);
    end
    figure(3); % Graphic for background subtraction.
    subplot(2, 2, 1:2); % Show difference.
    imagesc(VidRaw - VidQuad);
    title('Difference');
    set(gca, 'XTick', [], 'YTick', []); % Remove axes and labels.
    subplot(2, 2, 3);
    imagesc(VidRaw); % Show raw.
    title('Raw image');
    set(gca, 'XTick', [], 'YTick', []); % Remove axes and labels.
    subplot(2, 2, 4);
    imagesc(VidQuad); % Show background subtraction.
    title('Quadratic in y subtracted image');
    set(gca, 'XTick', [], 'YTick', []); % Remove axes and labels.
    set(gcf, 'Color', 'w');
    writeVideo(QuadVideo, getframe(gcf)); % Write to video.
    
    %% Gabor filtering.
    for j = 0:GStep:180 % Step through filter angles.
        % Perform Gabor filtering.
        [GMag, GPhase] = imgaborfilt(VidFrame(:, :, 1), GFilter, j);
        % Get total signal at particular angle as output.
        GArray((j / GStep) + 1, (i / VidStep) + 1) = sum(sum(GMag));
        figure(1);
        subplot(2, 2, 3);
        imshow(GMag, [], 'InitialMagnification', 'Fit'); % Show Gabor magnitude.
        title('Gabor magnitude');
        subplot(2, 2, 4);
        imshow(GPhase, [], 'InitialMagnification', 'Fit'); % Show Gabor phase.
        title('Gabor phase');
        set(gcf, 'Color', 'w');
        writeVideo(GaborVideo, getframe(gcf));
    end
    GArray(:, (i / VidStep) + 1) = GArray(:, (i / VidStep) + 1)...
        / max(GArray(:, 1)); % Normalise data to first frame.    
    subplot(2, 2, 2);
    imagesc(GArray);
    title('Angular intensity');
    xlabel('Time [frame number]');
	ylabel('Gabor filter angle');
    set(gca, 'YTick', [1, (90 / GStep) + 1, size(GArray, 1)],...
        'YTickLabels', [0, 90, 180], 'YDir', 'Reverse');
%         LBPFrame = extractLBPFeatures(mag(:,:,1), 'Upright', false,...
%             'NumNeighbors', 24);
%         if i == 1
%             LBPRef = LBPFrame;
%         end
%         LBPSurf(:, i) = (LBPFrame - LBPRef).^2;
%         subplot(2, 2, 2);
%         surf(LBPSurf, 'EdgeColor', 'None');
%         imagesc(LBPSurf);
%         title('Binary pattern');
end

%% Save videos.
close(QuadVideo);
close(GaborVideo);

%% Show surface of angular Gabor intensity.
figure(2);
surf(GArray, 'EdgeColor', 'None');
view(gca,[-30 45]);
set(gca, 'ZTick', [], 'YTick', [1, 10, 18],...
    'YTickLabels', [0, 90, 180], 'YDir', 'Reverse');
set(gcf, 'Color', 'w');
xlabel('Time [frame number]');
ylabel('Gabor filter angle');
zlabel('Normalised Gabor magnitude');
