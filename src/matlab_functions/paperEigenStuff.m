% Damp
%% Damped Processing
% processEigenModesDamped("in/2d_tile_20by20/100by20/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])
% processEigenModesDampedPara("in/2d_tile_20by20/40by40/", "out/2d_damped_eigenStuff/", [1,.75,.5,.25, 0.1, 0.01])
% processEigenModesDampedPara("in/2d_lattice/", "out/2d_damped_eigenStuff/", [1,.75,.5,.25, 0.1, 0.01,.001, 0])
% processEigenModesDampedPara("in/2d_lattice/", "out/2d_damped_eigenStuff/", [10, 7.5, 5, 2.5],"periodic", true)
% processEigenModesDampedPara("in/2d_tile_20by20/40by40/", "out/2d_damped_eigenStuff/staging/", [],"periodic", true)

%% Damped Mode Data
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat") 

%% Damped Mode Density PDF
% 40 by 40
plotDampedModeDensityPDF(outData, [.005, .01, .04,.16], [1])
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
plotData = filterData(outData, 'pressure', .1, 'damping', .001)
plotRealImagEigenValues(plotData);
slopeLine('loglog' ,2, [1,3], 2E-3, 'TextLocation', [1.4, 2E-3])
slopeLine('loglog' ,.75, [3,24], 2E-2, 'TextLocation', [7, 3E-2])

% High pressure high damping
plotData = filterData(outData, 'pressure', .1, 'damping', .5)
plotRealImagEigenValues(plotData);
slopeLine('loglog' ,2, [1,3], 6E-2, 'TextLocation', [1.4, 6E-2])
slopeLine('loglog' ,1, [3,24], 6E-1, 'TextLocation', [7, 3E-1])

% low pressure high damping
plotData = filterData(outData, 'pressure', .001, 'damping', .5)
plotRealImagEigenValues(plotData);
slopeLine('loglog' ,.25, [1,3], 6E-1, 'TextLocation', [1.4, 2E-1])
slopeLine('loglog' ,1, [3,24], 6E-1, 'TextLocation', [7, 3E-1])
%% damped sandbox
matMass = eye(2)
% matDamp = [1 -1; -1 1]
matDamp = zeros(2)
matSpring = [2 -1; -1 2]

[eigVec, eigVal] = eig(matSpring, matMass)
[polyVec, polyVal] = polyeig(matSpring, matDamp, matMass)


% atan2(imag(polyVec(1,4),real(polyVec(1,4))))


% Keeping one of the two conjugate pair eigenValues
keep = imag(polyVal) >0
shapes = polyVec(:,keep) % correspoding eigenVectors

for k = 1:size(shapes,2) % go through each column
    [~, idx]   = max(abs(shapes(:,k))) % get the largest eigenvector this column
    % normalize mode "k" by the largest, turns it into a unit phasor
    phase      = shapes(idx,k) / abs(shapes(idx,k))
    shapes(:,k)= shapes(:,k) / phase        % rotate all the eigenvectors by the same amount
end

% 3. take the real part
shapes = real(shapes)


modeVector = eigenVectors(:,1);
% This is what is in plotEigenmode
[~, idx]   = max(abs(modeVector)); % get the largest eigenvector this column
% Turn largest eigenvector into a unit phasor
unitPhasor =  modeVector(idx) / abs(modeVector(idx))
modeVector= modeVector(:) / unitPhasor % rotate all the eigenvectors by the same amount