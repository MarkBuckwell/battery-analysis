function [] = MyPlotOptions(FontSize)
% Function to quickly set my own default plot options.
%   Because Matlab won't let me switch the box off as default.
set(gcf, 'Color', 'White');
set(gca, 'FontSize', FontSize, 'FontName', 'Arial', 'LineWidth', 2,...
    'TickDir', 'Out', 'Box', 'Off', 'LineWidth', 2, 'Colormap', jet(64));
end

