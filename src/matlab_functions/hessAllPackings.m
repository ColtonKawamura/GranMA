function hessAllPackings(file_name_list, path_to_store)
    % hessAllPackings computes the eigenvectors and eigenvalues of the Hessian 
    % matrix for each packing in the specified list of files and saves them in a 
    % separate file.
    % 
    % Inputs:
    %   file_name_list: A cell array or string array of full paths to the *.mat files.
    %   path_to_store: Path where the eigenvectors and eigenvalues should be saved.

    % Loop through each file in the list
    for i = 1:length(file_name_list)
        % Get the current file name
        file_name = file_name_list{i};
        
        try
            % Try to load the required variables from the .mat file
            load(file_name, 'x', 'y', 'Dn', 'K', 'Ly', 'Lx');
            
            % Prepare positions and radii
            positions = [x', y'];
            radii = Dn ./ 2;
            
            % Clean the data using the cleanRats function (you need to define cleanRats)
            [positions, radii] = cleanRats(positions, radii, K, Ly, Lx);
            
            % Compute the Hessian matrix using hess2d (you need to define hess2d)
            Hessian = hess2d(positions, radii, K, Ly, Lx);
            
            % Calculate the eigenvalues and eigenvectors of the Hessian
            [eigen_vectors, eigen_values] = eig(Hessian);
            
            % Define the output filename (using the original file name without the .mat extension)
            [~, filename, ~] = fileparts(file_name); % Get the base file name
            output_file_name = fullfile(path_to_store, ['hess_', filename, '_eig.mat']); % Add "hess_" prefix
            
            % Save the eigenvectors and eigenvalues in the specified location using MAT-file version 7.3
            save(output_file_name, 'eigen_vectors', 'eigen_values', '-v7.3');
            
            % Clear all variables to free memory after processing this file
            clear x y Dn K Ly Lx positions radii Hessian eigen_vectors eigen_values;
            
            % Display progress
            fprintf('Processed file: %s\n', file_name);
            
        catch ME
            % If an error occurs (e.g., missing variables), print a warning and continue
            fprintf('Error processing file: %s. Skipping this file.\n', file_name);
            fprintf('Error message: %s\n', ME.message);
        end
    end
end
