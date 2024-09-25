function [slope1, slope2 ] = getSlopeK(x, y)
    % Remove NaN or invalid values from the data
    validData = ~isnan(y) & ~isnan(x); % Assuming NaN values might exist
    x = x(validData);
    y = y(validData);

    % Perform k-means clustering to separate the two lines
    k = 2;  % Number of clusters (two lines)
    data = [x(:), y(:)];
    [idx, ~] = kmeans(data, k);

    % Separate the data into two clusters
    x1 = x(idx == 1);
    y1 = y(idx == 1);
    x2 = x(idx == 2);
    y2 = y(idx == 2);

    % Fit a line to the first cluster using polyfit
    p1 = polyfit(x1, log(y1), 1);  % Linear fit on log-scale for y1
    slope1 = p1(1);

    % Fit a line to the second cluster using polyfit
    p2 = polyfit(x2, log(y2), 1);  % Linear fit on log-scale for y2
    slope2 = p2(1);

    % Plot the original data
    figure;
    scatter(x, y, 'o');
    hold on;

    % Plot the cleaned data
    scatter(x1, y1, 'filled');
    scatter(x2, y2, 'filled');

    % Plot the fitted lines
    x_fit = linspace(min(x), max(x), 100);
    plot(x_fit, exp(polyval(p1, x_fit)), 'r', 'LineWidth', 2, 'DisplayName', 'Fitted Line 1');
    plot(x_fit, exp(polyval(p2, x_fit)), 'g', 'LineWidth', 2, 'DisplayName', 'Fitted Line 2');

    % Add legend and labels
    xlabel('X-axis');
    ylabel('Y-axis (log scale)');
    legend('Original Data', 'Cluster 1', 'Cluster 2', 'Fitted Line 1', 'Fitted Line 2');
    set(gca, 'YScale', 'log');

    % Display slopes of the two lines
    fprintf('Slope of Line 1: %.4f\n', slope1);
    fprintf('Slope of Line 2: %.4f\n', slope2);

    hold off;
end
