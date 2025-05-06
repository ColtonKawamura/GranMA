% Damp
%% Damped Processing
processEigenModesDamped("in/2d_tile_20by20/100by20/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01, 0.001])
processEigenModesDampedPara("in/2d_tile_20by20/40by40/", "out/2d_damped_eigenStuff/", [1, 0.1, 0.01])


%% Damped Mode Density PDF
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1600_40by56_K100_M1.mat", "outData")
load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1891_100by28_K100_M1.mat", "outData")

% 40 by 40
plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,0, [.1,1], .45, 'TextLocation', [.5, .5])
slopeLine('loglog' ,1, [.1,1.5], .09, 'TextLocation', [.75, .1])

%  40 by 40 undamped test
slopeLine('loglog' ,1, [.025,.25], .62, 'TextLocation', [.08, 1])
slopeLine('loglog' ,4, [.1,.25], .15, 'TextLocation', [.2, .15])

% 100 by 28
plotDampedModeDensityPDF(outData, [.2, .01, .001], [.1])
slopeLine('loglog' ,1, [.01,.1], 2.5, 'TextLocation', [.03, 3.5])
slopeLine('loglog' ,4, [.031,.15], 1.3, 'TextLocation', [.075, 1])


%% Damped Eigen Vectors
    load("out/2d_damped_eigenStuff/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat", "outData"); 
    outData = orderPolyEig(outData);
    plotData = filterData(outData, 'pressure', .1, 'damping', .001)
    x = plotData.positions{1}(:, 1);
    y = plotData.positions{1}(:, 2);
    eigenVectors = plotData.eigenVectors{1};
    modeToPlot = 1;
    plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);

%% UnDamped Mode Density PDF
% 40 by 40
fileNameList = [
    "in/2d_tile_20by20/40by40/2D_N1600_P0.1_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.05_Width40_Seed1.mat",
    "in/2d_tile_20by20/40by40/2D_N1600_P0.01_Width40_Seed1.mat",
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
load("out/junkyard/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat") % test - damped algo with undamped data
% ---- Maybe point where needs own function ------
positions = [x',y']; 
radii = Dn/2;  
Hessian = hess2d(positions, radii', K, Ly, Lx);
[eigenVectors, eigenValues ] = eig(Hessian);
% -------------------------------------
modesToPlot = 1:10;
for i = 1:length(modesToPlot)
    modeToPlot = modesToPlot(i);
    plotEigenmode(x', y', eigenVectors, modeToPlot);
    input('Press Enter to continue...');
end


 % Gaussian
 sim2dGauss(100, 1, .0001, 1, 5000, 0.1, 10, 1, "in/2d_5wide_1000long/", "out/junk_yard")



%  Small Packing for Debugging
%  Process the small packing
% processEigenModesDampedPara("in/2d_damped_eigen_small/", "out/junkyard/", [0])
load("out/junkyard/2D_damped_eigenstuff_N14_5by4_K100_M1.mat", "outData"); % Small packing
outData = orderPolyEig(outData);
plotData = filterData(outData, 'pressure', .1, 'damping',0)
% [~, idx] = sort(abs(imag(plotData.eigenValues{1})));
% plotData.eigenVectors{1} = plotData.eigenVectors{1}(:, idx);
x = plotData.positions{1}(:, 1);
y = plotData.positions{1}(:, 2);
eigenVectors = plotData.eigenVectors{1};
modeToPlot = 1;
plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);

% This verifies that the spring matrix is correct
% It does what processEigenModesDamped does when you do matSpringDampMass
load("in/2d_damped_eigen_small/2D_N14_P0.1_Width3_Seed1.mat")
positions = [x',y']; 
radii = Dn'/2;  
[Hessian, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, 0,1 );
[eigenVectors, eigenValues] = polyeig(Hessian, matDamp, matMass);
% [eigenVectors, eigenValues ] = eig(Hessian); % this plots correct
% sort the eigenValues and eigenVectors
[~, idx] = sort(abs(imag(eigenValues)));
eigenVectors = eigenVectors(:, idx);

x = positions(:, 1);
y = positions(:, 2);
modeToPlot = 1;
modesToPlot = 1:size(eigenVectors,2);
for i = 1:length(modesToPlot)
    modeToPlot = modesToPlot(i);
    plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);
    input('Press Enter to continue...');
end
plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', true);

% Undamped
load("in/2d_damped_eigen_small/2D_N14_P0.1_Width3_Seed1.mat")
positions = [x',y']; 
radii = Dn'/2;  
Hessian = hess2d(positions, radii, K, Ly, Lx);
[eigenVectors, eigenValues ] = eig(Hessian);
x = positions(:, 1);
y = positions(:, 2);
modeToPlot = 1;
plotEigenmode(x, y, eigenVectors, modeToPlot, 'damped', false);

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