function nearestValue = findNear(list, value)
    % findNear: Finds the value in the list closest to the given scalar value
    % Inputs:
    %   list - A vector of numerical values
    %   value - A scalar value to compare
    % Output:
    %   nearestValue - The value in the list closest to the input value

    % Validate inputs
    if ~isvector(list)
        error('Input "list" must be a vector.');
    end
    if ~isscalar(value)
        error('Input "value" must be a scalar.');
    end

    % Calculate the absolute differences between each element and the value
    differences = abs(list - value);

    % Find the index of the smallest difference
    [~, index] = min(differences);

    % Retrieve the closest value
    nearestValue = list(index);
end
