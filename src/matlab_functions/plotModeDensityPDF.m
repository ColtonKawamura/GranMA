function plotModeDensityPDF(file_name_list, process_data)
% Plots the PDF of Eigen Modes. Input is list of packings like this:
% 
% file_name_list1 = [
%     "in/packings_processed_eig_PDF/2D_N6400_P0.1_Width40_Seed1.mat", 
%     "in/packings_processed_eig_PDF/2D_N6400_P0.01_Width40_Seed1.mat", 
%     "in/packings_processed_eig_PDF/2D_N6400_P0.001_Width40_Seed1.mat"
% ];
% 
% If the packings have not been pre-processed, set process_data = true

pressure_list = zeros(size(file_name_list));
for i = 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name);
    load(file_name, "P");
    pressure_list(i) = P;
end
figure
for i = 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name)
    positions = [x',y'];  % Assuming x, y are in the loaded .mat file
    radii = Dn./2;  % Assuming Dn is in the loaded .mat file

    if process_data == true
        Hessian = hess2d(positions, radii, K, Ly, Lx);
        [eigen_vectors, eigen_values ] = eig(Hessian);
    end

    [edges, normalized_counts] = modeDensity(eigen_values);
    [~, marker_color] = normVarColor(pressure_list, P, 1);

    plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', sprintf('$ P = %.3f $', P));
    xlabel('eigen frequencies (edges)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('counts', 'Interpreter', 'latex', 'FontSize', 20)
    legend('show', 'Interpreter', 'latex');
    title(sprintf('$L_x$ by $L_y$: %.2f by %.2f', Lx, Ly), 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, "XScale", "log")
    set(gca, "YScale", "log")
    grid on; 
    hold on;
end

end