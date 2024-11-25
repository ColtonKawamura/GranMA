% This is a script not a function.

%  Make the packing
N = 30;
Dn_list = 1+0.4*(rand([1,N])>0.5); 
packper103f(1,1,Dn_list,0,1);
% file_name = sprintf('in/2D_Repeating/2D_N%d_P0_Width0_Seed1.mat', N); % fully periodic

N = 5000;
% file_name = sprintf('in/2D_N5000_P0.1_Width5_Seed1.mat'); % regular packing
N = 26;
file_name = sprintf('in/2D_N26_P0.1_Width5_Seed1.mat'); % regular packing

load(file_name)

[Hessian, eigen_values, eigen_vectors] = HessYale(x, y, Dn, N, Ly, K, Lx);

% Plot Mode
mode_to_plot = 1; % should just be the column
figure;
plotEigenmode(x', y', eigen_vectors, mode_to_plot)
eigen_values(1:10, 1:10)