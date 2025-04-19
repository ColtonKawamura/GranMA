function filteredData = filterData3d(data, varargin)
    % Filters a struct array based on spatial limits within fields or closest scalar field values.
    % Uses a two-pass approach: calculates combined spatial masks first, then applies them.
    %
    % Args:
    %   data: A struct array (e.g., data loaded from .mat file).
    %   varargin: Pairs of arguments for filtering:
    %       - Spatial Filtering: { [lower, upper], 'y' } or { [lower, upper], 'z' }
    %         Filters vector fields within each struct element based on initial y or z coordinates.
    %       - Scalar Field Filtering: { value, 'field_name_string' }
    %         Filters the entire struct array, keeping elements where 'field_name'
    %         is closest to 'value'. Field must be numeric.
    %
    % Returns:
    %   filtered_data: The filtered struct array.


    if ~isscalar(data) || ~isstruct(data)
        error('Input data must be a scalar struct.');
    end

    if mod(length(varargin), 2) ~= 0
        error('Filtering arguments must be provided in field_name/value pairs.');
    end


    allFieldNames = fieldnames(data);

    for i = 1:2:length(varargin)
        fieldName = varargin{i}
        value = varargin{i+1}
        whos value
        % fullFieldData = data.(fieldName); % all the values for fieldName for all simulations

        % absDif = abs(fullFieldData - value); 
        % minDiff = min(absDif);
        % tolerance = minDef * 1e-9;
        % matchingIndex = find(absDif <= minDiff+ tolerance);
    end

end % End function