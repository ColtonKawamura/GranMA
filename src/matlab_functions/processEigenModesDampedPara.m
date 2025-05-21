function processEigenModesDampedPara(in_path, out_path, dampingConstants, options)

% Processes eigenmodes for damped granular packings with parallel computation
%
% This function loads granular packing data from .mat files, calculates the 
% eigenmodes with various damping constants, and saves the results. Uses 
% parallel processing to handle multiple files simultaneously.
%
% Parameters:
% in_path - Directory containing .mat files with granular packing data
% out_path - Directory to save the processed results
% dampingConstants - Array of damping constants to apply to each packing
% options.periodic - Boolean flag for periodic boundary conditions (default: false)

arguments
    in_path (1,1) string
    out_path (1,1) string
    dampingConstants (1,:) double
    options.periodic (1,1) logical = false
end
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
outData.Ly             = zeros(numTotalCombinations, 1);
outData.Lx             = zeros(numTotalCombinations, 1);
outData.radii              = cell(numTotalCombinations, 1);
outData.positions              = cell(numTotalCombinations, 1);
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
    radii = loadedVars.Dn' / 2;
    [positions, radii] = cleanRats(positions, radii, loadedVars.K, loadedVars.Ly, loadedVars.Lx);

    % Inner loop (serial)
    for j = 1:numDamping
        dampingConstant = dampingConstants(j);
        fprintf('CPU %d processing pressure %f with damping %d \n', i, loadedVars.P, dampingConstant);

        if options.periodic
            [Hessian, matDamp, matMass] = matSpringDampMass(positions, radii, loadedVars.Ly, loadedVars.Lx, dampingConstant, "periodic", true);
        else
            [Hessian, matDamp, matMass] = matSpringDampMass(positions, radii, loadedVars.Ly, loadedVars.Lx, dampingConstant);
        end

        [eigenVectors_j, eigenValues_j] = polyeig(Hessian, matDamp, matMass);

        % Store outData for this (i, j) pair in a temporary struct
        data_j = struct();
        data_j.pressure = loadedVars.P;
        data_j.damping = dampingConstant;
        data_j.eigenVectors = eigenVectors_j;
        data_j.eigenValues = eigenValues_j;
        data_j.radii = radii;
        data_j.Ly = loadedVars.Ly;
        data_j.Lx = loadedVars.Lx;
        data_j.positions = positions;

        % Put this struct into the cell for the current file i
        workerData_i{j} = data_j;
    end

    % Assign the cell containing outData for file i to the overall temporary cell
    allWorkerData{i} = workerData_i;
end

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
        outData.radii{idx} = data_j.radii;
        outData.Ly(idx) = data_j.Ly;
        outData.Lx(idx) = data_j.Lx;
        outData.positions{idx} = data_j.positions;
    end
end


% --- Saving ---
first_valid_idx = find(outData.pressure ~= 0 | ~cellfun('isempty', outData.eigenValues), 1, 'first');
if ~isempty(first_valid_idx)
    N_save = length(outData.positions{first_valid_idx});
    Lx_save = outData.Lx(first_valid_idx);
    Ly_save = outData.Ly(first_valid_idx);
    K_save = 100; % Placeholder - adjust as needed if K varies and needs to be stored/retrieved
    mass_save = 1;
    filename_output = sprintf("2D_damped_eigenstuff_N%d_%dby%d_K%d_M%d.mat", N_save, round(Lx_save), round(Ly_save), K_save, mass_save);
else
    error('No valid data found to save.');
end

save_path = fullfile(out_path, filename_output);
% Remove unneeded index fields before saving
outData = rmfield(outData, {'source_file_index', 'damping_index'});
save(save_path, 'outData', '-v7.3');
fprintf("Saved %d outData to: %s\n", numTotalCombinations, save_path);

end