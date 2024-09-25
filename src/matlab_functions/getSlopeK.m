function [slope1, slope2] = getKSlope(x, y)
    % Remove NaN or invalid values from the data
    validData = ~isnan(y) & ~isnan(x); % Assuming NaN values might exist
    x = x(validData);
    y = y(validData);

    % Log-transform the y-values to handle the log scale
    log_y = log(y);

    % Sort the y-values and calculate the differences
    [sorted_log_y, sort_idx] = sort(log_y);
    diffs = diff(sorted_log_y);

    % Find the largest gap in the sorted log(y) values
    [~, max_gap_idx] = max(diffs);

    % Use the mid-point of the largest gap as the threshold
    threshold_log_y = (sorted_log_y(max_gap_idx) + sorted_log_y(max_gap_idx + 1)) / 2;

    % Assign points to clusters based on the threshold
    idx = log_y > threshold_log_y; % Logical index for top cluster

    % Separate the data into two clusters (top and bottom)
    x1 = x(idx);      % Top cluster (higher y-values)
    y1 = y(idx);
    x2 = x(~idx);     % Bottom cluster (lower y-values)
    y2 = y(~idx);

    % Fit a line to the first cluster (top)
    p1 = polyfit(x1, log(y1), 1);  % Linear fit on log-scale for y1
    slope1 = p1(1);

    % Fit a line to the second cluster (bottom)
    p2 = polyfit(x2, log(y2), 1);  % Linear fit on log-scale for y2
    slope2 = p2(1);

    % Plot the original data
    figure;
    scatter(x, y, 'o', 'DisplayName', 'Original Data');
    hold on;

    % Plot the separated clusters
    scatter(x1, y1, 'filled', 'DisplayName', 'Top Cluster');
    scatter(x2, y2, 'filled', 'DisplayName', 'Bottom Cluster');

    % Plot the fitted lines
    x_fit = linspace(min(x), max(x), 100);
    plot(x_fit, exp(polyval(p1, x_fit)), 'r', 'LineWidth', 2, 'DisplayName', 'Fitted Line 1 (Top)');
    plot(x_fit, exp(polyval(p2, x_fit)), 'g', 'LineWidth', 2, 'DisplayName', 'Fitted Line 2 (Bottom)');

    % Add legend and labels
    xlabel('X-axis');
    ylabel('Y-axis (log scale)');
    legend('Location', 'best');
    set(gca, 'YScale', 'log');

    % Display slopes of the two lines
    fprintf('Slope of Line 1 (Top): %.4f\n', slope1);
    fprintf('Slope of Line 2 (Bottom): %.4f\n', slope2);

    hold off;
end
