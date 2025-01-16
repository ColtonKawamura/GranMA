function output_data = findInStruct(your_struct, searched_fieldname, searched_values, output_fieldname)
    % findInStruct - This function searches for values in a specified field of a structure 
    % and retrieves the corresponding data from another field.
    %
    % Inputs:
    %   your_struct - A structure array containing data.
    %   searched_fieldname - The name of the field to search in each element of the structure.
    %   searched_values - A vector or array of values to search for in searched_fieldname.
    %   output_fieldname - The name of the field from which to retrieve data when a match is found.
    %
    % Output:
    %   output_data - The data found in the output_fieldname for the matched elements.
    % Example:
    % findInStruct(test_struct, "combo", [.2, .02], "data")
    
    % Initialize an empty array to store the matching data
    output_data = [];
    
    % Loop through the structure array
    for idx = 1:length(your_struct)
        % Get the current value in the searched field
        current_value = your_struct(idx).(searched_fieldname);
        
        % Check if current_value matches the searched_values
        % Here we assume searched_values is a vector or array
        if isequal(current_value, searched_values)
            % If a match is found, retrieve the corresponding data from the output_fieldname
            if isfield(your_struct(idx), output_fieldname)
                % Store the data from the output_fieldname
                output_data = [output_data; your_struct(idx).(output_fieldname)];
            else
                % If the output_fieldname does not exist in the current struct, print a warning
                warning('Field "%s" not found in structure at index %d', output_fieldname, idx);
            end
        end
    end
end
