function slope = getSlope(x, y)

    fit_line = polyfit(x, y, 1);  
    y_fitted = polyval(fit_line, x);  

    residuals = abs(y - y_fitted);
    threshold = 2 * std(residuals);

    valid_idx = residuals < threshold;
    x_clean = x(valid_idx);  
    y_clean = y(valid_idx);

    fit_line_clean = polyfit(x_clean, y_clean, 1);
    slope = fit_line_clean(1);

    disp(['Slope of the line (after removing outliers): ', num2str(slope)]);

    figure;
    plot(x, y, 'o');
    hold on;
    plot(x_clean, y_clean, 'bo', 'MarkerFaceColor', 'b');  
    plot(x, polyval(fit_line_clean, x), 'r-', 'LineWidth', 2);  
    legend('Original Data', 'Cleaned Data', 'Fitted Line');
    xlabel('X-axis');
    ylabel('Y-axis');
    hold off;

end
