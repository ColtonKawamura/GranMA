function slope = getSlopeLog(x, y)

    fit_line = polyfit(x, log(y), 1);  % Logarithmic fit for semilogy plot
    slope = fit_line(1);
    
    
    figure;
    semilogy(x, y, 'o');
    hold on;
    semilogy(x, exp(polyval(fit_line, x)), 'r-', 'LineWidth', 2); % Plot the fitted line
    legend('Data', 'Fitted Line');
    xlabel('X-axis');
    ylabel('Y-axis (log scale)');
    hold off;
    
    end
    