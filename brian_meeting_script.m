% Scipt for meeting with Brian

% _______________________________________________________________________________
% Visualizing  Damped Eigenmodes

% Parameters
load("in/2d_damped_eigen_small/2D_N14_P0.1_Width3_Seed1.mat")
load("in/2d_eigen_mode_test/2D_N400_P0.1_Width10_Seed1.mat")
damping_constant = .1
mode_to_plot = 1;

% processing 
positions = [x', y'];
radii = Dn ./ 2;
[positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
mass = 1;
[matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
[eigen_vectors, eigen_values] = polyeig(matSpring, matDamp, matMass);
plotEigenmodeDamped(x, y, eigen_vectors, mode_to_plot)
figure
loglog(abs(imag(eigen_values)), -real(eigen_values), 'o')
xlabel('frequency', 'FontSize', 14)
ylabel('damping', 'FontSize', 14)
grid on

% _______________________________________________________________________________
%  Damped Density of Modes 40 long 12 wide

pressure_list = [.2, .01, .001];
damping_list = [.1]; % usually just pick one
plotDampedModeDensityPDF("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N400_K100_M1.mat", pressure_list, damping_list)

% _______________________________________________________________________________
%  Non-damped Density of Modes 160 and 45 wide
file_name_list = [
    "in/packings_processed_eig_PDF/2D_N6400_P0.1_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N6400_P0.01_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N6400_P0.001_Width40_Seed1.mat"
];

plotModeDensityPDF(file_name_list, false) % If the packings have not been pre-processed, set process_data = true