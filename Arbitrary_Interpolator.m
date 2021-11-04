function [Ix1, Iy1, Ix2, Iy2, interpOut] = Arbitrary_Interpolator(x1, y1, x2, y2)
%% This function will interpolate two datasets to equivalent sampling intervals.
% I.e. input two pairs of column vectors, [x1, y1] and [x2, y2], and the
% function will produce resampled versions whose datapoints may
% subsequently be compared directly at the same indices.

%% Sort out the starting workspace.
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

%% Resample each vector to increase resolution by nUp/nDen.
figure;
nUp = 1; % Upsample scalar.
nDen = 1; % Upsample denominator.
upFilter = 1; % Upsampling filter.
upB = 1; % Upsampling filter shape factor.
x1Up = resample(x1, nUp, nDen, upFilter, upB);
y1Up = resample(y1, nUp, nDen, upFilter, upB);
subplot(2, 2, 1);
hold on
plot(x1, y1);
plot(x1Up, y1Up);

x2Up = resample(x2, nUp, nDen, upFilter, upB);
y2Up = resample(y2, nUp, nDen, upFilter, upB);
subplot(2, 2, 2);
hold on
plot(x2, y2);
plot(x2Up, y2Up);

%% Resample each vector to evenly spaced increments.
x1Scale = resample(x1Up, x1Up);
y1Scale = resample(y1Up, x1Up);
subplot(2, 2, 3);
hold on
plot(x1, y1);
plot(x1Scale, y1Scale);

x2Scale = resample(x2Up, x2Up);
y2Scale = resample(y2Up, x2Up);
subplot(2, 2, 4);
hold on
plot(x2, y2);
plot(x2Scale, y2Scale);

%% Downsample longer dataset, interpolating to shorter datapoints.
if numel(x1Scale) > numel(x2Scale)
    Ix1 = interp1(x1Scale, x1Scale, x2Scale, 'Nearest');
    Iy1 = interp1(x1Scale, y1Scale, x2Scale, 'Nearest');
    Ix2 = x2Scale;
    Iy2 = y2Scale;
else
    Ix2 = interp1(x2Scale, x2Scale, x1Scale, 'Nearest');
    Iy2 = interp1(x2Scale, y2Scale, x1Scale, 'Nearest');
    Ix1 = x1Scale;
    Iy1 = y1Scale;
end

interpOut = [Ix1, Iy1, Ix2, Iy2];
open interpOut;
figure;
hold on
plot(x1, y1);
plot(Ix1, Iy1);
yyaxis right;
plot(x2, y2);
plot(Ix2, Iy2);