function output_data = findInStruct(your_struct, searched_fieldnames, searched_values, varargin)
    % findInStruct - This function searches for values in multiple specified fields of a structure 
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
    
    % Loop through the structure array to find matching entries
    for idx = 1:length(your_struct)
        % Initialize a flag to track if all search conditions are met
        match_found = true;
        
        % Check each input field and its corresponding value
        for i = 1:length(searched_fieldnames)
            current_field = searched_fieldnames{i};
            expected_value = searched_values{i};
            
            % Check if the current field matches the expected value
            if ~isequal(your_struct(idx).(current_field), expected_value)
                match_found = false;
                break;  % No need to check further fields, as one mismatch is enough
            end
        end
        
        % If all conditions are met, retrieve the corresponding data from the output fields
        if match_found
            for i = 1:length(varargin)
                output_fieldname = varargin{i};
                
                % Check if the field exists
                if isfield(your_struct(idx), output_fieldname)
                    % Store the data from the output field
                    output_data{i} = your_struct(idx).(output_fieldname);
                else
                    % If the output field is not found, give a warning
                    warning('Field "%s" not found in structure at index %d', output_fieldname, idx);
                end
            end
        end
    end
end
