function processEigenModesDamped(in_path, out_path, damping_constants)

% pressures = ["0.001", "0.01", "0.1"];
% out_path = "out/"
% damping_constants = [1, 0.1, 0.01, 0.001];
% processEigenModesDamped("in/2d_eigen_mode_test", "out/2d_damped_eigenStuff", [1, 0.1, 0.01, 0.001]);
% processEigenModesDamped("in/2d_damped_eigen_small", "out/junkyard", [1, 0.1, 0.01, 0.001]);

data = struct();
filenameList = dir(fullfile(in_path, '*.mat'));

for i = 1:length(filenameList)
    filename = fullfile(in_path, filenameList(i).name);

    try
        load(filename, 'x', 'y', 'Dn', 'K', 'Ly', 'Lx', 'P', 'N');
    catch
        warning("File %s does not contain the expected variables. Skipping...", filename);
        continue;
    end

    positions = [x', y'];
    radii = Dn ./ 2;
    [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
    mass = 1;

    for j = 1:length(damping_constants)
        damping_constant = damping_constants(j);
        
        [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
        
        [eigen_vectors, eigen_values] = polyeig(matSpring, matDamp, matMass);
        
        data.pressure(i,1) = P;
        data.damping(i,1) = damping_constant;
        data.eigen_vectors{i,1} = eigen_vectors;
        data.eigen_values{i,1} = eigen_values;
        data.diameter{i,1} = Dn;
        data.Ly(i,1) = Ly;
        data.Lx(i,1) = Lx;
        data.x{i,1} = x;
        data.y{i,1} = y;
    end
end

filename_output = string(sprintf("2D_damped_eigenstuff_N%d_%dby%d_K%d_M%d.mat", N, Lx, round(Ly), K, mass));
save_path = fullfile(out_path, filename_output);
save(save_path, 'data', '-v7.3'); % need this for file sizes larger than 2G's
display("Saved to: " + save_path);
data
end