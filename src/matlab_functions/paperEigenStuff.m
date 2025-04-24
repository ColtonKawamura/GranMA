% Damped
%% Processing
processEigenModesDamped("in/2d_tile_20by20/20by40/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])

%% Plotting
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1200_60by23_K100_M1.mat", "data")

plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,1.3, [.003,.1], 1.5, 'TextLocation', [.018, 3.5])

% UnDamped
file_name_list = [
    "in/packings_processed_eig_PDF/2D_N1600_P0.1_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N1600_P0.01_Width40_Seed1.mat", 
    "in/packings_processed_eig_PDF/2D_N1600_P0.001_Width40_Seed1.mat"
];

fileNameList = [
    "in/2d_tile_20by20/80by20/2D_N1200_P0.1_Width20_Seed1.mat",
    "in/2d_tile_20by20/80by20/2D_N1200_P0.01_Width20_Seed1.mat",
    "in/2d_tile_20by20/80by20/2D_N1200_P0.001_Width20_Seed1.mat"
];
 plotModeDensityPDF(fileNameList, true)

 % Gaussian
 