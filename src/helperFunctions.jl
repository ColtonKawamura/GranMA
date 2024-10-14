export
    fitLogLine,
    cleanVector

function fitLogLine(x, y)
    X = hcat(ones(length(x)), x)  # Create matrix for linear regression [1 x]
    coeffs = X \ log.(y)          # Perform linear regression on log(y)

    # Return coefficients [intercept, slope]
    return coeffs
end

function cleanVector(x, y)
    mean_y = mean(y)
    std_y = std(y)

    # Create a mask that keeps points within 2 standard deviations
    mask = abs.(y .- mean_y) .<= .5 * std_y

    # Apply the mask to both x and y to keep corresponding x values
    filtered_x = x[mask]
    filtered_y = y[mask]

    return filtered_x, filtered_y
end
