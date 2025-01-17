function plotDampedModeDensityPDF(struct_file, pressure_list, damping_list)
% Plots the PDF of Damped Eigen Modes. Input is a struct with fields :
% 
% pressure, damping, eigen_values, eigen_vectors
% 
% struct_file = *.mat file processed by ______ 
% pressure_list = vector of pressures you want to be plotted
% damping_list = damping list you want to have plotted
% 
% Note: usually you want multiple presssures and a single damping constant
% More note: this assumes that the eigen_values = damping + i frequency
% plotDampedModeDensityPDF("eigen_results.mat", [.1, .01, .001], [.01])

load(struct_file)
damping_constant = damping_list(1); % This is temporary unless we need a reason to plot multiple damping on one plot

figure
for i = 1:length(pressure_list)
    pressure = pressure_list(i);

    % positions = [x',y'];  % Assuming x, y are in the loaded .mat file
    % radii = Dn./2;  % Assuming Dn is in the loaded .mat file

    eigen_values = findInStruct(results, {'pressure', 'damping'}, {pressure, damping_constant}, 'eigen_values'); 
    eigen_values = eigen_values{1}; % need to do this because the struct saves it as an 1-deep array

    [edges, normalized_counts] = modeDensity(abs(imag(eigen_values)));
    [~, marker_color] = normVarColor(pressure_list, pressure, 1);

    % plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', sprintf('$ P = %.3f $', pressure));
    plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color, 'Color', marker_color);
    xlabel('eigen frequencies (edges)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('counts', 'Interpreter', 'latex', 'FontSize', 20)
    % legend('show', 'Interpreter', 'latex');
    % title(sprintf('$L_x$ by $L_y$: %.2f by %.2f', Lx, Ly), 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, "XScale", "log")
    set(gca, "YScale", "log")
    grid on; 
    hold on;
end

end