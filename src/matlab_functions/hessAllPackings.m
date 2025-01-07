function hessAllPackings(path_to_packings, path_to_store)
    % hessAllPackings computes the eigenvectors and eigenvalues of the Hessian 
    % matrix for each packing in the specified directory and saves them in a 
    % separate file.
    % 
    % Inputs:
    %   path_to_packings: Path where the input *.mat files are located.
    %   path_to_store: Path where the eigenvectors and eigenvalues should be saved.

    % Get list of all *.mat files in the specified directory
    file_list = dir(fullfile(path_to_packings, '*.mat'));
    
    % Loop through each file in the list
    for i = 1:length(file_list)
        % Construct the full file name
        file_name = fullfile(path_to_packings, file_list(i).name);
        
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
        output_file_name = fullfile(path_to_store, [filename, '_eig.mat']);
        
        % Save the eigenvectors and eigenvalues in the specified location
        save(output_file_name, 'eigen_vectors', 'eigen_values');
        clear x y Dn K Ly Lx positions radii Hessian eigen_vectors eigen_values; % If I get an OOM fault I'm going to throw a chair
        
        fprintf('Processed file: %s\n', file_name);
    end
end
