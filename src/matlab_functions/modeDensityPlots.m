file_name_list1 = [
    "in/2D_N1600_P0.1_Width40_Seed1.mat", 
    "in/2D_N1600_P0.01_Width40_Seed1.mat", 
    "in/2D_N1600_P0.001_Width40_Seed1.mat"
];

file_name_list2 = [
    "in/2D_N1600_P0.1_Width20_Seed1.mat", 
    "in/2D_N1600_P0.01_Width20_Seed1.mat", 
    "in/2D_N1600_P0.001_Width20_Seed1.mat"
];

file_name_list3 = [
    "in/2D_N900_P0.1_Width10_Seed1.mat", 
    "in/2D_N900_P0.01_Width10_Seed1.mat", 
    "in/2D_N900_P0.001_Width10_Seed1.mat"
];

file_name_list = file_name_list1
file_name_list = file_name_list2
file_name_list = file_name_list3
% file_name_list = [file_name_list1, file_name_list2, file_name_list3]; % concatenate as string array

% just get the pressure for normalization 
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
    % Ly=Ly*10;  % Optional adjustment if necessary
    positions = [x',y'];  % Assuming x, y are in the loaded .mat file
    radii = Dn./2;  % Assuming Dn is in the loaded .mat file
    % Skip the calculations for positions, radii, Hessian, and eigenvalues

    % Hessian = hess2d(positions, radii, K, Ly, Lx);
    % [eigen_vectors, eigen_values ] = eig(Hessian);
    [edges, normalized_counts] = modeDensity(eigen_values);
    [~, marker_color] = normVarColor(pressure_list, P, 1);

    plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', sprintf('$ P = %.3f $', P));
    xlabel('eigen frequencies (edges)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('counts', 'Interpreter', 'latex', 'FontSize', 20)
    legend('show', 'Interpreter', 'latex');
    title(sprintf('$L_x$ by $L_y$: %.2f by %.2f', Lx, Ly), 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, "XScale", "log")
    grid on; 
    hold on;
end



%__________________________________________________________

file_name_list1 = [
    "in/2D_N1600_P0.1_Width40_Seed1.mat", 
    "in/2D_N1600_P0.1_Width20_Seed1.mat", 
    "in/2D_N900_P0.1_Width10_Seed1.mat"
];

% P = 0.01 files
file_name_list2 = [
    "in/2D_N1600_P0.01_Width40_Seed1.mat", 
    "in/2D_N1600_P0.01_Width20_Seed1.mat", 
    "in/2D_N900_P0.01_Width10_Seed1.mat"
];

% P = 0.001 files
file_name_list3 = [
    "in/2D_N1600_P0.001_Width40_Seed1.mat", 
    "in/2D_N1600_P0.001_Width20_Seed1.mat", 
    "in/2D_N900_P0.001_Width10_Seed1.mat"
];
file_name_list = file_name_list1
file_name_list = file_name_list2
file_name_list = file_name_list3

% file_name_list = [file_name_list1, file_name_list2, file_name_list3]; % concatenate as string array

% just get the width for normalization 
width_list = zeros(size(file_name_list));
for i = 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name);
    load(file_name, "P");
    width_list(i) = Ly;
end
figure
for i = 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name)
    % Ly=Ly*10;  % Optional adjustment if necessary
    positions = [x',y'];  % Assuming x, y are in the loaded .mat file
    radii = Dn./2;  % Assuming Dn is in the loaded .mat file
    % Skip the calculations for positions, radii, Hessian, and eigenvalues

    % Hessian = hess2d(positions, radii, K, Ly, Lx);
    % [eigen_vectors, eigen_values ] = eig(Hessian);
    [edges, normalized_counts] = modeDensity(eigen_values);
    [~, marker_color] = normVarColor(width_list, Ly, 1);

    plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color, 'Color', marker_color, 'DisplayName', sprintf('$ Ly = %.3f $', Ly));
    xlabel('eigen frequencies (edges)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('counts', 'Interpreter', 'latex', 'FontSize', 20)
    legend('show', 'Interpreter', 'latex');
    title(sprintf('Pressure: %.3f ', P), 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, "XScale", "log")
    grid on; 
    hold on;
end


%---------------------------------------------------------
file_name_list1 = [
    "in/2D_N6400_P0.1_Width40_Seed1.mat", 
    "in/2D_N6400_P0.01_Width40_Seed1.mat", 
    "in/2D_N6400_P0.001_Width40_Seed1.mat"
];

file_name_list2 = [
    "in/2D_N3600_P0.1_Width20_Seed1.mat", 
    "in/2D_N3600_P0.01_Width20_Seed1.mat", 
    "in/2D_N3600_P0.001_Width20_Seed1.mat"
];

file_name_list3 = [
    "in/2D_N6400_P0.1_Width80_Seed1.mat", 
    "in/2D_N6400_P0.01_Width80_Seed1.mat", 
    "in/2D_N6400_P0.001_Width80_Seed1.mat"
];

file_name_list4= [
    "in/2D_N3600_P0.1_Width10_Seed1.mat", 
    "in/2D_N3600_P0.01_Width10_Seed1.mat", 
    "in/2D_N3600_P0.001_Width10_Seed1.mat"
];

file_name_list = file_name_list1
file_name_list = file_name_list2
file_name_list = file_name_list3
file_name_list = file_name_list4
% file_name_list = [file_name_list1, file_name_list2, file_name_list3]; % concatenate as string array

% just get the pressure for normalization 
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
    % Ly=Ly*10;  % Optional adjustment if necessary
    positions = [x',y'];  % Assuming x, y are in the loaded .mat file
    radii = Dn./2;  % Assuming Dn is in the loaded .mat file
    % Skip the calculations for positions, radii, Hessian, and eigenvalues

    % Hessian = hess2d(positions, radii, K, Ly, Lx);
    % [eigen_vectors, eigen_values ] = eig(Hessian);
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