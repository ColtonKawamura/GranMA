function slopeLine(type, slope, x_bounds, y_center)
% SLOPELINE Plot a line with a given slope and center point.
%   SLOPELINE(TYPE, SLOPE, X_BOUNDS, Y_CENTER) plots a line with the given
%   slope and center point. The line is plotted over the x values in
%   X_BOUNDS. The type of plot is determined by the string TYPE, which can
%   be 'linear', 'semilog', or 'loglog'. The center point is given by
%   Y_CENTER.

x = x_bounds;
hold on;
% Calculate y values based on the type of plot
switch type
    case 'linear'
        y = slope * (x - mean(x_bounds)) + y_center;
    case 'semilog'
        y = slope * (log10(x) - mean(log10(x_bounds))) + y_center;
    case 'loglog'
        y = 10.^(slope * (log10(x) - mean(log10(x_bounds)))) * y_center;
    otherwise
        error('Invalid type. Use "linear", "semilog", or "loglog".');
end

plot(x, y, 'k', 'LineWidth', 2, 'HandleVisibility', 'off'); % Exclude from legend
hold on;
text(mean(x_bounds), mean(y) - 0.05 * range(y), ['$' sprintf('%.2f', slope) '$'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 14, 'Interpreter', 'latex');

