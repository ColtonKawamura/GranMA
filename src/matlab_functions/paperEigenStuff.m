% Damped
%% Processing
% processEigenModesDamped("in/2d_tile_20by20/20by40/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])
processEigenModesDampedPara("in/2d_tile_20by20/tiles/40by40/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01])

%% Plotting
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1600_40by56_K100_M1.mat", "outData")

plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,1, [.01,.1], .1, 'TextLocation', [.03, 3.5])
slopeLine('loglog' ,2, [.01,1], .6, 'TextLocation', [.018, 3.5])

% UnDamped
fileNameList = [
    "in/2d_tile_20by20/40by40/2D_N1600_P0.1_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.01_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.001_Width40_Seed1.mat"
];
plotModeDensityPDF(fileNameList, true)
slopeLine('loglog' ,1, [.08,1.3], .08, 'TextLocation', [.5, .07])


 % Plot Eigen Vectors Damped
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat", "outData"); 
plotData = filterData(outData, 'pressure', .1, 'damping', .1)
x = plotData.positions{1}(:, 1);
y = plotData.positions{1}(:, 2);
eigenVectors = plotData.eigenVectors{1};
modeToPlot = 1;
plotEigenmodeDamped(x, y, eigenVectors, modeToPlot);

 % Gaussian
 sim2dGauss(100, 1, .0001, 1, 5000, 0.1, 10, 1, "in/2d_5wide_1000long/", "out/junk_yard")