function processEigenModesDampedPara(in_path, out_path, dampingConstants)

% pressures = ["0.001", "0.01", "0.1"];
% out_path = "out/"
% dampingConstants = [1, 0.1, 0.01, 0.001];
% processEigenModesDamped("in/2d_eigen_mode_test", "out/2d_damped_eigenStuff", [1, 0.1, 0.01, 0.001]);
% processEigenModesDamped("in/2d_damped_eigen_small", "out/junkyard", [1, 0.1, 0.01, 0.001]);

filenameList = dir(fullfile(in_path, '*.mat'));
numFiles = length(filenameList);
numDamping = length(dampingConstants);
numTotalCombinations = numFiles * numDamping;

% --- Pre-allocation for the final structure ---
if numFiles == 0
    warning('No .mat files found in %s', in_path);
    return;
end

outData = struct();
outData.pressure       = zeros(numTotalCombinations, 1);
outData.damping        = zeros(numTotalCombinations, 1);
outData.eigenVectors  = cell(numTotalCombinations, 1);
outData.eigenValues   = cell(numTotalCombinations, 1);
outData.diameter       = cell(numTotalCombinations, 1);
outData.Ly             = zeros(numTotalCombinations, 1);
outData.Lx             = zeros(numTotalCombinations, 1);
outData.x              = cell(numTotalCombinations, 1);
outData.y              = cell(numTotalCombinations, 1);
outData.source_file_index = zeros(numTotalCombinations, 1);
outData.damping_index   = zeros(numTotalCombinations, 1);

% --- Create a temporary cell array to store outData from workers ---
% Each cell will hold the outData for one file i (across all j)
allWorkerData = cell(numFiles, 1);

% --- Parallel Loop ---
parfor i = 1:numFiles
    filename = fullfile(in_path, filenameList(i).name);

    % Create a temporary structure for outData of this specific file i
    workerData_i = cell(numDamping, 1);

    loadedVars = load(filename, 'x', 'y', 'Dn', 'K', 'Ly', 'Lx', 'P', 'N');

    positions = [loadedVars.x', loadedVars.y'];
    radii = loadedVars.Dn ./ 2;
    [positions, radii] = cleanRats(positions, radii, loadedVars.K, loadedVars.Ly, loadedVars.Lx);
    mass = 1;

    % Inner loop (serial)
    for j = 1:numDamping
        dampingConstant = dampingConstants(j);
        fprintf('CPU %d processing pressure %f with damping %d \n', i, loadedVars.P, dampingConstant);

        [matSpring, matDamp, matMass] = matSpringDampMass(positions, radii, loadedVars.K, loadedVars.Ly, loadedVars.Lx, dampingConstant, mass);
        [eigenVectors_j, eigenValues_j] = polyeig(matSpring, matDamp, matMass);

        % Store outData for this (i, j) pair in a temporary struct
        data_j = struct();
        data_j.pressure = loadedVars.P;
        data_j.damping = dampingConstant;
        data_j.eigenVectors = eigenVectors_j;
        data_j.eigenValues = eigenValues_j;
        data_j.diameter = loadedVars.Dn;
        data_j.Ly = loadedVars.Ly;
        data_j.Lx = loadedVars.Lx;
        data_j.x = loadedVars.x;
        data_j.y = loadedVars.y;
        data_j.source_file_index = i;
        data_j.damping_index = j;

        % Put this struct into the cell for the current file i
        workerData_i{j} = data_j;
    end

    % Assign the cell containing outData for file i to the overall temporary cell
    allWorkerData{i} = workerData_i;
end

% --- Assemble outData after parfor loop ---
% --- Assemble outData after parfor loop ---
for i = 1:numFiles
    workerData_i = allWorkerData{i}; % Get the cell array of results for file i

    for j = 1:numDamping
        data_j = workerData_i{j}; % Get the struct for damping j

        % Calculate the linear index for the pre-allocated outData structure
        idx = (i - 1) * numDamping + j;

        % Assign data from the temporary struct to the final pre-allocated struct
        outData.pressure(idx) = data_j.pressure;
        outData.damping(idx) = data_j.damping;
        outData.eigenVectors{idx} = data_j.eigenVectors;
        outData.eigenValues{idx} = data_j.eigenValues;
        outData.diameter{idx} = data_j.diameter;
        outData.Ly(idx) = data_j.Ly;
        outData.Lx(idx) = data_j.Lx;
        outData.x{idx} = data_j.x;
        outData.y{idx} = data_j.y;
        outData.source_file_index(idx) = data_j.source_file_index; % This should be i
        outData.damping_index(idx) = data_j.damping_index;     % This should be j
    end
end


% --- Saving ---
first_valid_idx = find(outData.pressure ~= 0 | ~cellfun('isempty', outData.eigenValues), 1, 'first');
if ~isempty(first_valid_idx)
    N_save = length(outData.x{first_valid_idx});
    Lx_save = outData.Lx(first_valid_idx);
    Ly_save = outData.Ly(first_valid_idx);
    K_save = 100; % Placeholder - adjust as needed if K varies and needs to be stored/retrieved
    mass_save = 1;
    filename_output = sprintf("2D_damped_eigenstuff_N%d_%dby%d_K%d_M%d_Combined.mat", N_save, Lx_save, round(Ly_save), K_save, mass_save);
else
    filename_output = '2D_damped_eigenstuff_NoResults.mat';
end

save_path = fullfile(out_path, filename_output);
save(save_path, 'outData', '-v7.3');
fprintf("Saved %d outData to: %s\n", numTotalCombinations, save_path);

end