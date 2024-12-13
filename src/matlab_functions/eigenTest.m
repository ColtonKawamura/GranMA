% This is a script not a function.

file_name = sprintf('in/2D_N5000_P0.01_Width5_Seed1.mat'); % regular packing
% file_name = sprintf('in/2D_N5000_P0.1_Width5_Seed1.mat'); % regular packing
% file_name = sprintf('in/2D_N5000_P0.046416_Width5_Seed2.mat')

load(file_name)
positions = [x',y'];
radii = Dn./2;
[positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
Hessian = hess2d(positions, radii, K, Ly, Lx);
% Hessian = TESTME(positions, radii, K, Ly, Lx);
[eigen_vectors, eigen_values ] =  eig(Hessian);
 % Plot Mode
% mode_to_plot = 100; % should just be the column
% figure;
% plotEigenmode(x', y', eigen_vectors, mode_to_plot)
% eigen_values(1:10, 1:10)

modes_to_plot = [1:10:200];
for i = 1:length(modes_to_plot)
    mode_to_plot = modes_to_plot(i);
    plotEigenmode(x', y', eigen_vectors, mode_to_plot)
    pause
end


%% Density of States
bin_size = 1;
bins = min(eigen_values(:)):bin_size:max(eigen_values(:));
counts = histcounts(eigen_values, bins);
histogram(sqrt(eigen_values), 'BinWidth', bin_size, 'Normalization', 'pdf')
set(gca,'YScale','log')

% Highest mode
[~, mode_bin_index] = max(counts);
mode_bin_center = (bins(mode_bin_index) + bins(mode_bin_index + 1)) / 2 % the one in the middle 

%% Debeye Plot
eigen_values_diag = diag(eigen_values);  % eigen values ar ethe diags

sqrt_eigen_values = sqrt(eigen_values_diag);

sorted_sqrt_eigen_values = sort(sqrt_eigen_values);

i = (1:length(sorted_sqrt_eigen_values))';  
N = length(positions);
figure;
semilogy(i./N, sorted_sqrt_eigen_values, 'o');
xlabel('$\frac{i}{N}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\log(\sqrt{\lambda})$', 'Interpreter', 'latex', 'FontSize', 20);
grid on;

%% Density of modes
file_name1 = "in/2D_N5000_P0.1_Width5_Seed1.mat"; % packing
file_name2 = "in/2D_N5000_P0.046416_Width5_Seed1.mat"; % regular packing
file_name3 = "in/2D_N5000_P0.021544_Width5_Seed1.mat"; % regular packing
file_name4 = "in/2D_N5000_P0.01_Width5_Seed1.mat"; % regular packing
file_name5 = "in/2D_N5000_P0.0046416_Width5_Seed1.mat"; % regular packing
file_name6 = "in/2D_N5000_P0.0021544_Width5_Seed1.mat"; % regular packing
file_name7 = "in/2D_N5000_P0.001_Width5_Seed1.mat"; % regular packing
file_name_list = [file_name1, file_name2, file_name3, file_name4, file_name5, file_name6, file_name7]; % concatenate as string array
    % just get the presssure for normalziation 
pressure_list = zeros(size(file_name_list));
for i = 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name);
    load(file_name, "P");
    pressure_list(i) = P;
end
for i= 1:length(file_name_list)
    file_name = file_name_list(i);
    load(file_name)
    positions = [x',y'];
    radii = Dn./2;
    [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
    Hessian = hess2d(positions, radii, K, Ly, Lx);
    [eigen_vectors, eigen_values ] =  eig(Hessian);
    [edges, normalized_counts] = modeDensity(eigen_values);
    [~, marker_color] = normVarColor(pressure_list, P, 1);
    plot(edges, normalized_counts, '-o', 'MarkerFaceColor', marker_color, 'DisplayName', sprintf('$ P = %.3f $', P)); 
    xlabel('eigen_frequencies (edges)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('counts', 'Interpreter', 'latex', 'FontSize', 20)
    legend('show', 'Interpreter', 'latex');
    grid on; 
    hold on;
end