function process_outputs_2d()
% Define the directory containing the .mat files
directory = 'outputs/seed_job_ellipse_edit/';

% Get a list of all .mat files in the directory and their metadata
mat_files = dir(fullfile(directory, '*.mat'));

% Initialize a table to store all data
data_table = table();
ellipse_data = {};

% Loop through each file and load data into the table
for k = 1:length(mat_files)
    % Construct the full file name and load the data
    file_name = fullfile(mat_files(k).folder, mat_files(k).name);
    file_data = load(file_name);
    
    if ~isempty(file_data.attenuation_fit_line_x)
        % Create a new row with the data
        new_row = table(file_data.pressure_dimensionless, ...
                        file_data.driving_angular_frequency_dimensionless, ...
                        file_data.gamma_dimensionless, ...
                        -file_data.attenuation_x_dimensionless, ...
                        file_data.attenuation_y_dimensionless, ...
                        -file_data.wavespeed_x, ...
                        -file_data.wavenumber_x_dimensionless, ...
                        file_data.seed, ...
                        file_data.mean_aspect_ratio, ...
                        file_data.mean_rotation_angles, ...
                        file_data.initial_distance_from_oscillation_output_x_fft(end), ...
                        file_data.input_pressure,...
                        'VariableNames', {'pressure', 'omega', 'gamma', 'attenuation_x', 'attenuation_y', 'wavespeed_x', 'wavenumber_x', 'seed', 'mean_aspect_ratio', 'mean_rotation_angles', 'attenuation_limit_x', 'input_pressure'});
                    
        % Append the new row to the table
        data_table = [data_table; new_row];
        
        
    end
end

% Add new columns to the table for plotting
data_table.omegagamma = data_table.omega .* data_table.gamma;
% This cuts off the highest 20 percent of the frequencies
% data_table = data_table(data_table.OmegaGamma < prctile(data_table.OmegaGamma, 90), :);
data_table = data_table(data_table.omegagamma < prctile(data_table.omegagamma, 100), :);
data_table.alphaoveromega = data_table.attenuation_x ./ data_table.omega;

writetable(data_table, 'crunched_data/2d_K100_ellipse_edits.csv')