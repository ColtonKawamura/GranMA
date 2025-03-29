function slopeLine(type, slope, x_bounds, y_center)
x = x_bounds;

% Calculate y values based on the type of plot
switch type
    case 'linear'
        y = slope * (x - mean(x_bounds)) + y_center;
    case 'semilog'
        y = slope * log10(x) + y_center;
    case 'loglog'
        y = 10.^(slope * log10(x) + y_center);
    otherwise
        error('Invalid type. Use "linear", "semilog", or "loglog".');
end

% Plot the line
plot(x, y, 'LineWidth', 2);
hold on;

