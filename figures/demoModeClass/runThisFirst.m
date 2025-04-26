load("2D_N400_P0.1_Width20_Seed1.mat"); 
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