% Damp
%% Damped Processing
% processEigenModesDamped("in/2d_tile_20by20/100by20/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])
% processEigenModesDampedPara("in/2d_tile_20by20/40by40/", "out/2d_damped_eigenStuff/", [1,.75,.5,.25, 0.1, 0.01])
% processEigenModesDampedPara("in/2d_lattice/", "out/2d_damped_eigenStuff/", [1,.75,.5,.25, 0.1, 0.01,.001, 0])
% processEigenModesDampedPara("in/2d_lattice/", "out/2d_damped_eigenStuff/", [10, 7.5, 5, 2.5],"periodic", true)
% processEigenModesDampedPara("in/2d_tile_20by20/40by40/", "out/2d_damped_eigenStuff/staging/", [],"periodic", true)

%% Damped Mode Data
load("out/2d_damped_eigenStuff/combined/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat") 

%% Damped Mode Density PDF
% 40 by 40
plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.25])
slopeLine('loglog' ,0, [.1,1], .45, 'TextLocation', [.5, .5])
slopeLine('loglog' ,1, [.1,1.5], .09, 'TextLocation', [.75, .1])
slopeLine('loglog' ,1/4, [5,1E3], 1, 'TextLocation', [100, .6]) % for the collapse

%% Damped Eigen Vectors
    load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat", "outData"); 
    outData = orderPolyEig(outData);
    plotData = filterData(outData, 'pressure', .1, 'damping', 2)
    x = plotData.positions{1}(:, 1);
    y = plotData.positions{1}(:, 2);
    eigenVectors = plotData.eigenVectors{1};
    modeToPlot = 1;
    plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);
    % for i = 1:10
    %     modeToPlot = i;
    %     plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);
    %     input('Press Enter to continue...');
    % end
    
%% UnDamped Mode Density PDF
% 40 by 40
fileNameList = [
    "in/2d_tile_20by20/40by40/2D_N1600_P0.1_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.075_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.05_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.025_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.01_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.005_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.001_Width40_Seed1.mat"
];
plotModeDensityPDF(fileNameList, true)
slopeLine('loglog' ,1, [.8,13], .008, 'TextLocation', [5, .007])
slopeLine('loglog' ,0, [1,10], .04, 'TextLocation', [4, .05])

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
% load("out/junkyard/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat") % test - damped algo with undamped data
% ---- Maybe point where needs own function ------
positions = [x',y']; 
radii = Dn/2;  
Hessian = hess2d(positions, radii', K, Ly, Lx);
[eigenVectors, eigenValues ] = eig(Hessian);
modesToPlot = 1:10;
for i = 1:length(modesToPlot)
    modeToPlot = modesToPlot(i);
    plotEigenmode(x', y', eigenVectors, modeToPlot);
    input('Press Enter to continue...');
end

%% Damped Imaginary vs Real Eigenvalues
outData = orderPolyEig(outData);

% High pressure low damping
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.005], "loglog", true);
slopeLine('linear' ,(.1*1.4)^-1, [1E-13,1E1], 40, 'TextLocation', [1, 2])
slopeLine('linear' ,(.01*1.4)^-1, [.01,1], 40, 'TextLocation', [1.4E-3, 2])
slopeLine('linear' ,(.001*1.4)^-1, [.01,.1], 40, 'TextLocation', [1.4E-3, 2])

% no scaling
slopeLine('linear' ,(.001*1.4)^-1, [0,5E-3], 1.5, 'TextLocation', [1.4E-3, 2])
slopeLine('linear' ,(.01*1.4)^-1, [0,5E-2], 1.5, 'TextLocation', [1.4E-3, 2])
slopeLine('linear' ,(.1*1.4)^-1, [0,.5], 1.5, 'TextLocation', [1.4E-3, 2])

slopeLine('loglog' ,2, [1E-2,1.1], 1E-5, 'TextLocation', [1.4E-3, 2])
slopeLine('loglog' ,1, [3E-1,1.3], 1E-3, 'TextLocation', [1.4E-3, 2])
slopeLine('loglog' ,.75, [3,24], 2E-2, 'TextLocation', [7, 3E-2])

% Plots for Brian paper

% Low Damping No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.001] );
slopeLine('linear' ,1.5, [0,2.5], 1.5, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.001], "loglog", true);
slopeLine('loglog' , .25, [7E-2,4E-1], 1.3, 'TextLocation', [1, 2])
slopeLine('loglog' , 1, [5E-1,5], 2, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [7E-2,4E-1], 1E-1, 'TextLocation', [1, 2])

% Low Damping Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.001], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.001], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])

% MedLow Damping No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.01] );
slopeLine('linear' ,1.5, [0,2.5], 1.5, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.01], "loglog", true);
slopeLine('loglog' , 1, [5E-1,5], 2, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,4E-1], 1E-1, 'TextLocation', [1, 2])

% MedLow Damping Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.01], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.01], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])

% Med Damping No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.05] );
slopeLine('linear' ,1.5, [0,2.5], 1.5, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.05], "loglog", true);
slopeLine('loglog' , 1, [5E-1,5], 2, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [7E-2,4E-1], 1E-1, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [7E-2,4E-1], 1.5, 'TextLocation', [1, 2])

% Med Damping Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.05], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.05], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])

% High Damping No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.1] );
slopeLine('linear' ,1.5, [0,2.5], 1.5, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.1], "loglog", true);
slopeLine('loglog' , 1, [5E-1,5], 2, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [7E-2,4E-1], 1E-1, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [7E-2,4E-1], 1.8, 'TextLocation', [1, 2])

% High Damping Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.1], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.1], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])

% OverDamped No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.25] );
slopeLine('linear' ,1.5, [0,2.5], 1.5, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.25], "loglog", true);
slopeLine('loglog' , 2, [1E-1,2], 1E-1, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [7E-2,4E-1], 1.8, 'TextLocation', [1, 2])

% Overdamped Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.25], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [.25], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])

% Way OverDamped No Scal
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [2] );
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [2], "loglog", true);

% Way Overdamped Scaled
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [2], 'scaling', true );
slopeLine('linear' ,1.5, [0,60], 30, 'TextLocation', [1.4E-3, 2])
plotRealImagEigenValues(outData, [.001, .005, .01, .05, .1, .5], [2], "loglog", true, 'scaling', true);
slopeLine('loglog' , 1, [2,80], 15, 'TextLocation', [1, 2])
slopeLine('loglog' , .25, [1,9], 30, 'TextLocation', [1, 2])
slopeLine('loglog' , 2, [1E-1,1], 1E-1, 'TextLocation', [1, 2])
%% SandBox

% Low Damping
plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.001])
slopeLine('loglog' ,0, [.1,1], .45, 'TextLocation', [.5, .5])
slopeLine('loglog' ,1, [.1,1.5], .09, 'TextLocation', [.75, .1])

plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.001], 'scaling', true)
% slopeLine('loglog' ,1/4, [5,1E3], 1, 'TextLocation', [100, .6]) % for the collapse

% Med Damping
plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.1])
slopeLine('loglog' ,0, [.1,1], .45, 'TextLocation', [.5, .5])
slopeLine('loglog' ,1, [.1,1.5], .09, 'TextLocation', [.75, .1])

plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.1], 'scaling', true)

% Over Damping
plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.25])
slopeLine('loglog' ,0, [.1,1], .45, 'TextLocation', [.5, .5])
slopeLine('loglog' ,1, [.1,1.5], .09, 'TextLocation', [.75, .1])
slopeLine('loglog' ,1/4, [5,1E3], 1, 'TextLocation', [100, .6]) % for the collapse

plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [.25], 'scaling', true)