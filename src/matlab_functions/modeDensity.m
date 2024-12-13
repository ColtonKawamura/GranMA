function [edges, normalized_counts] = modeDensity(eigen_values)
    eigen_values_diag = diag(eigen_values);  % eigen values ar ethe diags
    sqrt_eigen_values = sqrt(eigen_values_diag);
    sorted_sqrt_eigen_values = sort(sqrt_eigen_values);
    [counts, edges] = hist(sorted_sqrt_eigen_values./10, 200);
    norm = sum(counts)*(edges(2)-edges(1));
    normalized_counts = counts./norm;
    % figure; 
    % plot(edges, counts./norm, "-o"); grid on
end
