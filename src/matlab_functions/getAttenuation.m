function slope = getAttenuation(x, y)

    % Fit a single line to all the data in log scale
    fit_line = polyfit(x, log(y), 1);  % Logarithmic fit for semilogy plot
    slope = fit_line(1);  % Slope of the fitted line
    
    % Display the slope of the fitted line
    disp(['Slope of the line: ', num2str(slope)]);
    
    % Plot the original data and the fitted line
    figure;
    semilogy(x, y, 'o');  % Plot the original data
    hold on;
    semilogy(x, exp(polyval(fit_line, x)), 'r-', 'LineWidth', 2); % Plot the fitted line
    legend('Data', 'Fitted Line');
    xlabel('X-axis');
    ylabel('Y-axis (log scale)');
    hold off;
    
    end
    