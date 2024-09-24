function slope = getSlopeLog(x, y)

    fit_line = polyfit(x, log(y), 1);  
    y_fitted = exp(polyval(fit_line, x));  

    residuals = abs(y - y_fitted);
    threshold = 2 * std(residuals);

    valid_idx = residuals < threshold;
    x_clean = x(valid_idx);  
    y_clean = y(valid_idx);

    fit_line_clean = polyfit(x_clean, log(y_clean), 1);
    slope = fit_line_clean(1);

    figure;
    semilogy(x, y, 'o');
    hold on;
    semilogy(x_clean, y_clean, 'bo', 'MarkerFaceColor', 'b');  
    semilogy(x, exp(polyval(fit_line_clean, x)), 'r-', 'LineWidth', 2);  
    legend('Original Data', 'Cleaned Data', 'Fitted Line');
    xlabel('X-axis');
    ylabel('Y-axis (log scale)');
    hold off;

end
