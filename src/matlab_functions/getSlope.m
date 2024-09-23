function slope = getSlope(x, y)

fit_line = polyfit(x, y, 1);  % Simple linear fit
slope = fit_line(1); 

disp(['Slope of the line: ', num2str(slope)]);

figure;
plot(x, y, 'o');  % Plot the original data
hold on;
plot(x, polyval(fit_line, x), 'r-', 'LineWidth', 2); % Plot the fitted line
legend('Data', 'Fitted Line');
xlabel('X-axis');
ylabel('Y-axis');
hold off;

end
