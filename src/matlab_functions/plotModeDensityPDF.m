function plotModeDensityPDF(file_name_list, process_data)
% Plots the PDF of Eigen Modes. Input is list of packings like this:
% 
% file_name_list = [
%     "in/packings_processed_eig_PDF/2D_N6400_P0.1_Width40_Seed1.mat", 
%     "in/packings_processed_eig_PDF/2D_N6400_P0.01_Width40_Seed1.mat", 
%     "in/packings_processed_eig_PDF/2D_N6400_P0.001_Width40_Seed1.mat"
% ];
% 
% If the packings have not been pre-processed, set process_data = true
% Example: plotModeDensityPDF(file_name_list, false)


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
    positions = [x',y'];  
    radii = Dn./2;  

    if process_data == true
        Hessian = hess2d(positions, radii, K, Ly, Lx);
        [eigen_vectors, eigen_values ] = eig(Hessian);
    end

    [edges, normalized_counts] = modeDensity(eigen_values);
    % --- Verification Step ---
    bin_widths = diff(edges); % Calculate the width of each bin (Length N)
    num_bins = length(bin_widths);
    num_counts = length(normalized_counts);
    counts_for_bins = normalized_counts(1:num_bins); 
    integral_approx = sum(counts_for_bins .* bin_widths);
    fprintf('Area Under PDF = %.6f\n', integral_approx);
    % --- End Verification Step ---
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