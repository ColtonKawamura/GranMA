% Damp
%% Damped Processing
% processEigenModesDamped("in/2d_tile_20by20/100by20/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])
processEigenModesDampedPara("in/2d_tile_20by20/tiles/40by40/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01])

%% Damped Mode Density PDF
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1600_40by56_K100_M1.mat", "outData")
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1891_100by28_K100_M1.mat", "outData")

% 40 by 40
plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,1, [.01,.1], 2.5, 'TextLocation', [.03, 3.5])
slopeLine('loglog' ,4, [.031,.15], 1.3, 'TextLocation', [.075, 1])

% 100 by 28
plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,1, [.01,.1], 2.5, 'TextLocation', [.03, 3.5])
slopeLine('loglog' ,4, [.031,.15], 1.3, 'TextLocation', [.075, 1])

%% Damped Eigen Vectors
    load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat", "outData"); 
    load("out/junkyard/2D_damped_eigenstuff_N14_4.666667e+00by4_K100_M1.mat", "outData"); % Small packing
    plotData = filterData(outData, 'pressure', .1, 'damping', .1)
    x = plotData.positions{1}(:, 1);
    y = plotData.positions{1}(:, 2);
    eigenVectors = plotData.eigenVectors{1};
    modeToPlot = 1;
    plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);

%% UnDamped Mode Density PDF
% 40 by 40
fileNameList = [
    "in/2d_tile_20by20/40by40/2D_N1600_P0.1_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.01_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.001_Width40_Seed1.mat"
];
plotModeDensityPDF(fileNameList, true)
slopeLine('loglog' ,1, [.08,1.3], .08, 'TextLocation', [.5, .07])
slopeLine('loglog' ,0, [.1,1], .4, 'TextLocation', [.4, .5])

% 100 by 20
fileNameList = [
    "in/2d_tile_20by20/100by20/2D_N2000_P0.1_Width20_Seed1.mat",
    "in/2d_tile_20by20/100by20/2D_N2000_P0.01_Width20_Seed1.mat",
    "in/2d_tile_20by20/100by20/2D_N2000_P0.001_Width20_Seed1.mat"
];
plotModeDensityPDF(fileNameList, true)
slopeLine('loglog' ,1.2, [.08,1.3], .08, 'TextLocation', [.5, .07])
slopeLine('loglog' ,0, [.1,1], .4, 'TextLocation', [.4, .5])

%% Undamped - Plot Eigen Vectors
load("in/2d_tile_20by20/40by40/2D_N1600_P0.1_Width40_Seed1.mat"); 
load("in/2d_damped_eigen_small/2D_N14_P0.1_Width3_Seed1.mat")
% ---- Maybe point where needs own function ------
positions = [x',y']; 
radii = Dn/2;  
Hessian = hess2d(positions, radii', K, Ly, Lx);
[eigenVectors, eigenValues ] = eig(Hessian);
% -------------------------------------
modesToPlot = [1, 2, 3, 4];
for i = 1:length(modesToPlot)
    modeToPlot = modesToPlot(i);
    plotEigenmode(x', y', eigenVectors, modeToPlot);
    input('Press Enter to continue...');
end


 % Gaussian
 sim2dGauss(100, 1, .0001, 1, 5000, 0.1, 10, 1, "in/2d_5wide_1000long/", "out/junk_yard")