function [binned_semi_minor_values, binned_semi_major_values, binned_rot_angle_values, bin_centers] = process_ellipse_data(ellipse_stats_nonzero)
% Purpose - averages ellipse parameters semi-major, semi-minor axis, and rotation angle over a slices of volume of the simulation space in slices away from oscillating wall
%
% Format:   [binned_semi_minor_values, binned_semi_major_values, binned_rot_angle_values, bin_centers] = ...
%            process_ellipse(ellipse_stats_nonzero)
%
% Input:    ellipse_stats_nonzero          - semi-major, semi-minor axis, rotation angle, each particle's initial distance from oscillating wall
%
% Output:   binned_semi_minor_values       - average distance from center of particle time-space trajectory along ellipse minor axis over bin
%           binned_semi_major_values       - average distance from center of particle time-space trajectory along ellipse major axis over bin
%           binned_rot_angle_values        - average rotation angle of particle's ellipse relative to unit wave vector over bin
%           bin_centers                    - position of bins relative to oscillating wall
%
% Note:     This assumes each unit in the distance is a particle diameter

% % Find where the last "good" oscillation was based on the attenuation fft function and use that to shave off data further out
% index_position_limit = find(ellipse_stats_nonzero(:,4) == position_last_oscillation);
% ellipse_stats_nonzero = ellipse_stats_nonzero(1:index_position_limit,:)

% Extract distance from wall
distance_from_wall = ellipse_stats_nonzero(:, 4);

% Define the range of distances from the oscillating wall
min_distance = min(distance_from_wall);
max_distance = max(distance_from_wall);
distance_range = max_distance - min_distance;

% Define the number of "slices" to divide the particles up into bins, based on the distance range, rounding up
% This assumes each unit in the distance is a particle diameter
particle_diameter = 1;
num_bins = ceil(distance_range / particle_diameter);

% Define the bin edges
bin_edges = linspace(min_distance, max_distance, num_bins + 1);

% Initialize arrays to store the results
binned_data = cell(num_bins, 1);
bin_centers = zeros(num_bins, 1); % Used for plotting as the main parameter to mark each bin

% Group the data based on the distance
for current_bin = 1:num_bins
    % Find the indices of the data points that fall within the current bin edges
    bin_indices = distance_from_wall >= bin_edges(current_bin) & distance_from_wall < bin_edges(current_bin + 1);
    
    % Extract the data for the particles that fall in the current bin
    binned_data{current_bin} = ellipse_stats_nonzero(bin_indices, :);
    
    % Calculate the center of the bin (mean distance)
    bin_centers(current_bin) = mean(distance_from_wall(bin_indices));
end

% Initialize the ellipse parameters for fitting
binned_semi_major_values = zeros(num_bins, 1);
binned_semi_minor_values = zeros(num_bins, 1);
binned_rot_angle_values = zeros(num_bins, 1);

% Average the ellipse parameters over the bin for each particle
for current_bin = 1:num_bins
    if ~isempty(binned_data{current_bin})
        % Extract the elliptical parameters for the current bin
        binned_semi_major_values(current_bin) = mean(binned_data{current_bin}(:, 1)); % Semi-major axis mean across the bin
        binned_semi_minor_values(current_bin) = mean(binned_data{current_bin}(:, 2)); % Semi-minor axis mean across the bin
        binned_rot_angle_values(current_bin) = mean(abs(binned_data{current_bin}(:, 3))); % Rotation angle mean across the bin
    end
end

