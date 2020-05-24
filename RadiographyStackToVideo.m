clear
clc
addpath(cd);

%% Imports and converts a stack of images into a video.
% Images should be named in the same order that they should appear in the video.
% No compression applied, ignores image type and size, and converts to 8-bit.

[FileGroup, DataPath] = uigetfile('*.*', 'DialogTitle', 'Select image files:',...
    'MultiSelect', 'on' ) ; % Gets file names and location.
cd(DataPath); % Set active folder to curernt image directory.
NFC = length ( FileGroup ) ;  % Number of files chosen.

%% Reading and image conversion parameters.

FirstImage = 1; % Index of image to be first frame of video.
LastImage = NFC; % Index of image to be final frame of video.
ImageSkip = 1; % To skip mages between frames. Set to 1 to skip no images.
TimeIncrement = 0.0005 ; % Time between frames in seconds.
TimeBoxColour = 'White'; % Box colour for time stamp.
VideoRotation = 0; % Image to video rotation angle.
Precision = 1; % Image to double conversion precision.

% Start video writer.
RadiographyVideo = VideoWriter('RadiographyVideo.avi', 'Grayscale AVI');
% RadiographyVideo.LosslessCompression = 'True'; % No compression.
open(RadiographyVideo); % Open the video.

%% Read in the image stack.
warning('Off'); % Turn off imread warning on image type.
for i = FirstImage : ImageSkip : LastImage % Iterate from first to last chosen images.
    % Read, rotate and invert image.
    ReadImage = imcomplement(imrotate(imread(FileGroup{i}, VideoRotation));
    
    % Insert a text box with the frame time.
    FrameTime = (TimeIncrement * (i - 1));
	TimeString = sprintf('Time (s) = %.4f', FrameTime);
    TimePosition = [10 10]; % Position of time index box.
    ReadImage = insertText(ReadImage, TimePosition, TimeString,...
        'FontSize', 12, 'BoxColor', TimeBoxColour, 'BoxOpacity', 0.6);
    
    % Convert to uint8 and then to double, and just use R values from RGB.
    ProcessedImage = im2double(im2uint8((ReadImage(:,:,1))));
    % Append new frame to video.
    writeVideo(RadiographyVideo, ProcessedImage);
end
warning('On'); % Swtich warnings back on.
close(RadiographyVideo); % Save video.
