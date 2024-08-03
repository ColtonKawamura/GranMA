function [mean_aspect_ratio, mean_rotation_angles] =plot_ellipse_statistics(binned_semi_minor_values, binned_semi_major_values, binned_rot_angle_values, bin_centers, max_distance_from_wall)

valid_indices = bin_centers <= max_distance_from_wall;
mean_aspect_ratio = mean(binned_semi_minor_values(valid_indices) ./ binned_semi_major_values(valid_indices));
mean_rotation_angles = mean(abs(binned_rot_angle_values(valid_indices)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aspect and angle plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the fitted parameters as a function of distance
figure;
subplot(2, 1, 1);
plot(bin_centers, binned_semi_minor_values./binned_semi_major_values, '-o', 'DisplayName', ['$\frac{b}{a}$']);
ylabel('$\frac{b}{a}$', 'Interpreter', 'latex', 'FontSize', 20);
set(get(gca,'ylabel'),'rotation',0)
title(['$\frac{b}{a}_\mathrm{avg} =$' num2str(mean_aspect_ratio)], "Interpreter", "latex")
grid on;
xlim([0 max_distance_from_wall])
ylim([0 .8])

subplot(2, 1, 2)
plot(bin_centers, abs(binned_rot_angle_values), '-o', 'DisplayName', ['$\theta$']);
ylabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 20);
set(get(gca,'ylabel'),'rotation',0)
% legend('show', 'Interpreter', 'latex');
xlabel('Distance from Oscillating Wall (particle diameters)','Interpreter', 'latex', 'FontSize', 15);
grid on;
xlim([0 max_distance_from_wall])
ylim([0 30])
title(['$\theta_\mathrm{avg} =$' num2str(mean_rotation_angles)], "Interpreter", "latex")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log semi/major axis plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
semilogy(bin_centers, binned_semi_major_values, '-o', 'DisplayName', ['a']);
xlabel('Distance from Oscillating Wall (particle diameters)','Interpreter', 'latex', 'FontSize', 15);
ylabel('Displacement (particle diameters)', 'Interpreter', 'latex', 'FontSize', 15);

grid on;
hold on;
semilogy(bin_centers, binned_semi_minor_values, '-o', 'DisplayName', ['b']);
legend('show', 'Interpreter', 'latex');
xlim([0 max_distance_from_wall])