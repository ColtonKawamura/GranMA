function [normalized_variable, marker_color] = normVarColor(vector_list, value, use_log)
    % Normalize the input vector list using either a log or regular transformation,
    % and scale it between 0 and 1.
    %
    % Arguments:
    %   vector_list : the list of values to normalize
    %   value : the specific value to find the index of for coloring
    %   use_log : boolean (true for log normalization, false for regular normalization)
    %
    % Outputs:
    %   normalized_variable : the normalized vector
    %   marker_color : RGB color corresponding to the normalized value
    
    % Check if the input is a valid vector
    if isempty(vector_list) || numel(vector_list) < 2
        error('Input must be a non-empty vector with at least two elements.');
    end
    
    % Select normalization method
    if use_log
        % Apply the log transformation
        transformed_vector = log(vector_list);
    else
        % Use regular normalization (min-max scaling)
        transformed_vector = vector_list;
    end
    
    % Find the minimum and maximum of the transformed data
    min_transformed = min(transformed_vector);
    max_transformed = max(transformed_vector);
    
    % Normalize the transformed data between 0 and 1
    normalized_variable = (transformed_vector - min_transformed) / (max_transformed - min_transformed);
    
    % Find the index of the specific value in the vector list
    idx = find(vector_list == value, 1);
    
    % Check if the value is found
    if isempty(idx)
        error('The specified value is not in the vector list.');
    end
    
    % Get the normalized value for the given index
    normalized_value = normalized_variable(idx);
    
    % Assign a color based on the normalized value (RGB format)
    marker_color = [normalized_value, 0, 1 - normalized_value];
end
