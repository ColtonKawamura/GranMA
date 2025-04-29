function mean_stdevs = meanPhaseDev(x_values, y_values, bin_width)
    % Ensure x_values and y_values are sorted by x_values
    [x_values_sorted, sorted_indices] = sort(x_values);
    y_values_sorted = y_values(sorted_indices);
    
    % Make the bins
    min_x = min(x_values_sorted);
    max_x = max(x_values_sorted);
    bin_edges = min_x:bin_width:max_x;

    % Ensure the last edge includes max_x if it falls exactly on a bin edge
    if bin_edges(end) < max_x
        bin_edges = [bin_edges, max_x];
    end
    bin_stddevs = [];
    
    % Iterate over each bin to calculate standard deviation of y_values within the bin
    for i = 1:length(bin_edges) - 1
        bin_start = bin_edges(i);
        bin_end = bin_edges(i + 1);
        
        % Find indices of x_values within the current bin range
        % Include the start edge, exclude the end edge, except for the last bin
        if i < length(bin_edges) - 1
            bin_indices = find(x_values_sorted >= bin_start & x_values_sorted < bin_end);
        else
            % Include the end edge for the last bin
            bin_indices = find(x_values_sorted >= bin_start & x_values_sorted <= bin_end);
        end
        
        if ~isempty(bin_indices)
            bin_stddevs(end+1) = std(y_values_sorted(bin_indices));
        end
    end
    
    mean_stdevs = mean(bin_stddevs(~isnan(bin_stddevs)));
    
    % If all bins were empty or resulted in NaN, set as NaN
    if isempty(mean_stdevs)
        mean_stdevs = NaN;
end