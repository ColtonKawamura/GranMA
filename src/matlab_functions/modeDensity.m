function [edges, normalized_counts] = modeDensity(eigenValues, options)
    
    arguments
        eigenValues (:,:) double
        options.damped (1,1) logical = false
        options.logBins (1,1) logical = false
    end

    [s1, s2] = size(eigenValues);
    if s1 == 1 || s2 == 1
        eigenValues = eigenValues; % this is for the struct (damped flow)
    else
        eigenValues = diag(eigenValues); % this is for the case where the eigenValues are a matrix
    end

    if ~options.damped
        eigenValues = eigenValues( eigenValues > 0 & isfinite(eigenValues) ); % Need to go back and verify this
        eigenValues = sqrt(eigenValues);
    end

    sortedEigenValues = sort(eigenValues);
    

    % [counts, edges] = hist(sortedEigenValues./10, 30); % 100 bins, 10 is the natural frequency
    lowLimit = min(sortedEigenValues);
    highLimit = max(sortedEigenValues);
    nBins = 30;
    if options.logBins
        edges = logspace(log10(lowLimit), log10(highLimit), nBins+1); % nEdges gives nEdges-1 bins
    else
        edges = linspace(lowLimit, highLimit, nBins+1);
    end
    counts = histcounts(sortedEigenValues, edges);
    % norm = sum(counts)*(edges(2)-edges(1));
    binWidths = diff(edges);
    normalized_counts = counts./(sum(counts) * binWidths); % normalize the counts by the area under the curve
end
