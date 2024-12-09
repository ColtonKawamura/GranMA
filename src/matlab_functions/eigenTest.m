% This is a script not a function.

file_name = sprintf('in/2D_N1000_P0.1_Width4_Seed1.mat'); % regular packing

load(file_name)

% [Hessian, eigen_values, eigen_vectors] = HessYale(x, y, Dn, N, Ly, K, Lx);

positions = [x',y'];
Hessian = hess2d(positions, Dn/2, K, Ly);
[eigen_vectors, eigen_values ] =  eig(Hessian);

% Plot Mode
mode_to_plot = 100; % should just be the column
figure;
plotEigenmode(x', y', eigen_vectors, mode_to_plot)
eigen_values(1:10, 1:10)

modes_to_plot = [1,2,3,4,10, 20, 30, 50, 100];
for i = 1:length(modes_to_plot)
    mode_to_plot = modes_to_plot(i);
    plotEigenmode(x', y', eigen_vectors, mode_to_plot)
end


%% Density of States
bin_size = 1;
bins = min(eigen_values(:)):bin_size:max(eigen_values(:));
counts = histcounts(eigen_values, bins);
histogram(eigen_values, 'BinWidth', bin_size, 'Normalization', 'pdf')
set(gca,'YScale','log')

% Highest mode
[~, mode_bin_index] = max(counts);
mode_bin_center = (bins(mode_bin_index) + bins(mode_bin_index + 1)) / 2 % the one in the middle 