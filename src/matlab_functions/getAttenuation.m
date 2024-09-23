function [slope_lower, slope_upper] = getAttenuation(x, y)

% Find the median or mean value to use as a threshold to separate the two lines
threshold = mean(y);  % You can also use 'median(y)' if the data is skewed

% Separate the data into two sets: one above the threshold and one below
upper_line_idx = y > threshold;  % Logical index for upper line
lower_line_idx = y <= threshold; % Logical index for lower line

% Extract the x and y data for the upper line
x_upper = x(upper_line_idx);
y_upper = y(upper_line_idx);

% Extract the x and y data for the lower line
x_lower = x(lower_line_idx);
y_lower = y(lower_line_idx);

% Initialize the slopes to NaN
slope_upper = NaN;
slope_lower = NaN;

% Check if there are data points for the upper line
if ~isempty(x_upper) && ~isempty(y_upper)
    % Fit a linear polynomial to the upper line (first degree)
    p_upper = polyfit(x_upper, log(y_upper), 1);  % Using log(y_upper) for semilogy plot
    slope_upper = p_upper(1);  % Slope of the upper line
end

% Check if there are data points for the lower line
if ~isempty(x_lower) && ~isempty(y_lower)
    % Fit a linear polynomial to the lower line (first degree)
    p_lower = polyfit(x_lower, log(y_lower), 1);  % Using log(y_lower) for semilogy plot
    slope_lower = p_lower(1);  % Slope of the lower line
end

% Display the slopes
disp(['Slope of the upper line: ', num2str(slope_upper)]);
disp(['Slope of the lower line: ', num2str(slope_lower)]);

% Optionally, plot the fits on the original data for visualization
figure;
semilogy(x, y, 'o');  % Plot the original data
hold on;

% Plot the upper line fit if it exists
if ~isnan(slope_upper)
    semilogy(x_upper, exp(polyval(p_upper, x_upper)), 'r-', 'LineWidth', 2); % Upper line fit
end

% Plot the lower line fit if it exists
if ~isnan(slope_lower)
    semilogy(x_lower, exp(polyval(p_lower, x_lower)), 'g-', 'LineWidth', 2); % Lower line fit
end

legend('Data', 'Upper Line Fit', 'Lower Line Fit');
xlabel('X-axis');
ylabel('Y-axis (log scale)');
hold off;

end
