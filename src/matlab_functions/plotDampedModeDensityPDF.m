function plotDampedModeDensityPDF(struct_file, pressure_list, damping_list)
% Plots the PDF of Damped Eigen Modes. Input is a struct with fields :
% 
% pressure, damping, eigen_values, eigen_vectors
% 
% struct_file = *.mat file processed by processEigenModesDamped("in/2d_damped_eigen_small", "out/junkyard", [1, 0.1, 0.01, 0.001]); 
% pressure_list = vector of pressures you want to be plotted
% damping_list = damping list you want to have plotted
% 
% Note: usually you want multiple presssures and a single damping constant
% More note: this assumes that the eigen_values = damping + i frequency
% plotDampedModeDensityPDF("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N400_K100_M1.mat", [.2, .01, .001], [.1])
load(struct_file)

figure
for i = 1:length(pressure_list)
    pressure = pressure_list(i);

    for j = 1:length(damping_list)
        % positions = [x',y'];  % Assuming x, y are in the loaded .mat file
        % radii = Dn./2;  % Assuming Dn is in the loaded .mat file
        damping_constant = damping_list(j);
        eigen_values = findInStruct(results, {'pressure', 'damping'}, {pressure, damping_constant}, 'eigen_values'); 

        eigen_values = eigen_values{1}; % need to do this because the struct saves it as an 1-deep array

        [edges, normalized_counts] = modeDensity(abs(imag(eigen_values)));
        if length(pressure_list) == 1
            [~, marker_color] = normVarColor(damping_list, damping_constant, 1);
        else
            [~, marker_color] = normVarColor(pressure_list, pressure, 1);
        end

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

end