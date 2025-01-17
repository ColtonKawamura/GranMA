function processEigenModesDamped(in_path, out_path, damping_constants)

% pressures = ["0.001", "0.01", "0.1"];
% out_path = "out/"
% damping_constants = [1, 0.1, 0.01, 0.001];
% processEigenModesDamped("in/2d_eigen_mode_test", "out/2d_damped_eigenStuff", [1, 0.1, 0.01, 0.001]);

results = struct();
loop_counter = 0;
filename_list = dir(fullfile(in_path, '*.mat'));

for i = 1:length(filename_list)
    filename = fullfile(in_path, filename_list(i).name);
    load(filename)
    positions = [x', y'];
    radii = Dn ./ 2;
    [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
    mass = 1;

    for j = 1:length(damping_constants)
        loop_counter = loop_counter +1;
        damping_constant = damping_constants(j);
        
        [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, K, Ly, Lx, damping_constant, mass);
        
        [eigen_vectors, eigen_values] = polyeig(matSpring, matDamp, matMass);
        
        results(loop_counter).pressure = P;
        results(loop_counter).damping = damping_constant;
        results(loop_counter).eigen_vectors = eigen_vectors;
        results(loop_counter).eigen_values = eigen_values;
    end
end

% Save the results to a .mat file
% save('eigen_results.mat', 'results');
filename_output = string(sprintf("2D_damped_eigenstuff_N%d_K%d_M%d.mat", N, K, mass));
save_path = fullfile(out_path, filename_output);
save(save_path, 'results');

end