function newData = combineEigenData(dataOne, dataTwo)

    if isfield(dataOne, 'source_file_index')
        dataOne = rmfield(dataOne, {'source_file_index', 'damping_index'});
    end

    if isfield(dataTwo, 'source_file_index')
        dataTwo = rmfield(dataTwo, {'source_file_index', 'damping_index'});
    end

    % Combine the fields of dataOne and dataTwo
    fieldNames = fieldnames(dataOne);
    newData = struct();

    for i = 1:length(fieldNames)
        field = fieldNames{i};
        newData.(field) = [dataOne.(field); dataTwo.(field)];
    end
end