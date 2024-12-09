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

modes_to_plot = [1,2,3,4, 500];
for i = 1:length(modes_to_plot)
    mode_to_plot = modes_to_plot(i)
    plotEigenmode(x', y', eigen_vectors, mode_to_plot)
end
