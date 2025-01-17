function output_data = findInStruct(your_struct, searched_fieldnames, searched_values, varargin)
    % findInStruct - This function searches for the closest match based on multiple fields in a structure 
    % and retrieves the corresponding data from one or more output fields.
    %
    % Inputs:
    %   your_struct - A structure array containing data.
    %   searched_fieldnames - A cell array of field names to search in each element of the structure.
    %   searched_values - A cell array of values to search for in each corresponding searched_fieldname.
    %   varargin - A variable number of output field names to retrieve from the structure.
    %
    % Output:
    %   output_data - A cell array containing the data found in the output fields for the matched elements.
    %
    % Example:
    %   output_data = findInStruct(results, {'field2', 'field3'}, {0.2, 0.3}, 'field1')

    % Initialize output_data as a cell array to store data from multiple output fields
    output_data = cell(1, length(varargin));  
    
    % Initialize a variable to track the minimum total difference
    min_diff = Inf;
    best_match_idx = NaN;
    
    % Loop through the structure array to find the closest match
    for idx = 1:length(your_struct)
        % Initialize total difference for this element
        total_diff = 0;
        
        % Check each input field and its corresponding value
        for i = 1:length(searched_fieldnames)
            current_field = searched_fieldnames{i};
            expected_value = searched_values{i};
            
            % Calculate absolute difference for numeric fields
            if isnumeric(your_struct(idx).(current_field)) && isnumeric(expected_value)
                diff = abs(your_struct(idx).(current_field) - expected_value);
                total_diff = total_diff + diff;
            elseif ischar(your_struct(idx).(current_field)) && ischar(expected_value)
                % For strings, use strcmp (0 for no match, 1 for exact match)
                diff = double(~strcmp(your_struct(idx).(current_field), expected_value));
                total_diff = total_diff + diff;
            else
                % For other types, you can add more comparison methods
                if ~isequal(your_struct(idx).(current_field), expected_value)
                    total_diff = total_diff + 1;  % Arbitrary large penalty for mismatch
                end
            end
        end
        
        % Update the best match if the current total difference is smaller
        if total_diff < min_diff
            min_diff = total_diff;
            best_match_idx = idx;
        end
    end
    
    % If a closest match is found, retrieve the data for the output fields
    if ~isnan(best_match_idx)
        for i = 1:length(varargin)
            output_fieldname = varargin{i};
            
            % Check if the field exists
            if isfield(your_struct(best_match_idx), output_fieldname)
                % Store the data from the output field
                output_data{i} = your_struct(best_match_idx).(output_fieldname);
            else
                % If the output field is not found, give a warning
                warning('Field "%s" not found in structure at index %d', output_fieldname, best_match_idx);
            end
        end
    else
        warning('No match found.');
    end
end
