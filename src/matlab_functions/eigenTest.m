% This is a script not a function.

file_name = sprintf('in/2D_N5000_P0.01_Width5_Seed1.mat'); % regular packing
% file_name = sprintf('in/2D_N5000_P0.1_Width5_Seed1.mat'); % regular packing
% file_name = sprintf('in/2D_N5000_P0.046416_Width5_Seed2.mat')

load(file_name)
positions = [x',y'];
radii = Dn./2;
[positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
% [Hessian, Zn] = hess2d(positions, Dn/2, K, Ly, Lx);
Hessian = TESTME(positions, radii, K, Ly, Lx);
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
file_name = sprintf('in/2D_N5000_P0.1_Width5_Seed1.mat'); % regular packing
file_name = sprintf('in/2D_N5000_P0.01_Width5_Seed1.mat'); % regular packing
file_name = sprintf('in/2D_N5000_P0.001_Width5_Seed1.mat'); % regular packing
for i= 1:length(pressure_list)
    pressure = pressure_list(i);
    load(file_name)
    positions = [x',y'];
    radii = Dn./2;
    [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
    % [Hessian, Zn] = hess2d(positions, Dn/2, K, Ly, Lx);
    Hessian = TESTME(positions, radii, K, Ly, Lx);
    [eigen_vectors, eigen_values ] =  eig(Hessian);
    [edges, normalized_counts] = modeDensity(eigen_values)
    plot(edges, normalized_counts, "-o"); grid on; hold on
end



[counts, edges] = hist(sorted_sqrt_eigen_values./10, 200)
norm = sum(counts)*(edges(2)-edges(1))
figure; plot(edges, counts./norm, "-o"); grid on
% Next step is to do these plots for different presusres all normalized.