function [edges, normalized_counts] = modeDensity(eigen_values)
    [s1, s2] = size(eigen_values);
    if s1 == 1 || s2 == 1
        eigen_values_diag = eigen_values; % this is for the struct (damped flow)
    else
        eigen_values_diag = diag(eigen_values); % this is for the case where the eigen_values are a matrix
    end
    % eigen_values_diag = diag(eigen_values);  % eigen values ar ethe diags
    sqrt_eigen_values = sqrt(eigen_values_diag);
    sorted_sqrt_eigen_values = sort(sqrt_eigen_values);
    [counts, edges] = hist(sorted_sqrt_eigen_values./10, 100); % 100 bins, 10 is the natural frequency
    norm = sum(counts)*(edges(2)-edges(1));
    normalized_counts = counts./norm;
    % figure; 
    % plot(edges, counts./norm, "-o"); grid on
end
