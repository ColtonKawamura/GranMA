export
    fitLogLine,
    cleanVector,
    meanDistNeighbor,
    wrappedDistance,
    meanDeltaYNeighbor,
    meanPhaseDev,
    getMeanField

function getMeanField(filtered_data; plot = true)
    attenuation = filtered_data[1].alphaoveromega_x # will need to figure out how to adapt this for y later
    distance_from_wall = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    omega = filtered_data[1].omega # but this is dimensionelss

    A = filtered_data[1].pressure/100
    mean_field = A*exp.(-attenuation*omega*distance_from_wall)
    x_parra = filtered_data[1].initial_distance_from_oscillation_output_x_fft
    y = filtered_data[1].amplitude_vector_x    
    prime_field = abs.(y - mean_field)


    if plot==true
        mat"""
        figure
        scatter($(x_parra), $(prime_field), "o", "DisplayName", "x prime field")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end


    if plot==true
        mat"""
        scatter($(x_parra), $(y), "v", "DisplayName", "raw x-amplitude")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end
    
    if plot==true
        mat"""
        scatter($(x_parra), $(mean_field), "DisplayName", "mean field")
        hold on
        set(gca, 'YScale', 'log')
        grid on
        xlabel("\$ x \$", "Interpreter", 'latex', "FontSize", 15)
        ylabel("\$ A(x) \$", "Interpreter", 'latex', "FontSize", 15)
        set(get(gca, 'ylabel'), 'rotation', 0);
        box on
        hold on 
        legend("show")
        """
    end
    
    x_perp = filtered_data[1].initial_distance_from_oscillation_output_y_fft
    y = filtered_data[1].amplitude_vector_y
    if plot == true
        mat"""
        scatter($(x_perp), $(y), "*", "DisplayName", "raw y-amplitude")
        set(gca, 'YScale', 'log')
        grid on
        legend('show', 'Location', 'northeast', 'Interpreter', 'latex');
        set(gca, 'FontSize', 15);
        legend('FontSize', 15)
        """
    end
    return mean_field
end

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

function meanPhaseDev(x_values, y_values, bin_width)
   
    # Ensure x_values and y_values are sorted by x_values
    sorted_indices = sortperm(x_values)
    x_values = x_values[sorted_indices]
    y_values = y_values[sorted_indices]
    
    # Make the bins
    min_x = minimum(x_values)
    max_x = maximum(x_values)
    bin_edges = min_x:bin_width:max_x
    bin_stddevs = Float64[]
    
    # Iterate over each bin to calculate standard deviation of y_values within the bin
    for i in 1:length(bin_edges) - 1
        bin_start = bin_edges[i]
        bin_end = bin_edges[i + 1]
        
        # Find indices of x_values within the current bin range
        bin_indices = findall(x -> bin_start <= x < bin_end, x_values)
        
        if !isempty(bin_indices)
            push!(bin_stddevs, std(y_values[bin_indices]))
        end
    end
    mean_stdevs= mean(filter(!isnan, bin_stddevs))
    return mean_stdevs
end

function wrappedDistance(y1, y2)
    direct_dist = abs(y1 - y2)
    wrapped_dist = 2π - direct_dist
    return min(direct_dist, wrapped_dist)
end


function meanDeltaYNeighbor(x_values, y_values)
    n = length(x_values)
    closest_y_distances = zeros(n)
    
    for i in 1:n
        min_dist_y = 2π  # Initialize with a high value
        
        for j in 1:n
            if i != j
                dist_y = wrappedDistance(y_values[i], y_values[j])
                
                # Update the minimum distance in y direction if it's the smallest found
                if dist_y < min_dist_y
                    min_dist_y = dist_y
                end
            end
        end
        
        closest_y_distances[i] = min_dist_y
    end
    println(mean(closest_y_distances)) 
    return mean(closest_y_distances)
end

function meanDistNeighbor(x_values, y_values)
    points = hcat(x_values, y_values)
    n = size(points, 1)
    distances = zeros(n)
    
    for i in 1:n
        dist_to_others = zeros(n)
        
        for j in 1:n
            if i != j
                # Calculate Euclidean distance for x and wrapped distance for y
                dist_x = abs(x_values[i] - x_values[j])
                dist_y = wrappedDistance(y_values[i], y_values[j])
                dist_to_others[j] = sqrt(dist_x^2 + dist_y^2)
                # dist_to_others[j] = dist_y
            else
                dist_to_others[j] = Inf  # Exclude the point itself
                # dist_to_others[j] = NaN  # Exclude the point itself
            end
        end
        
        # Find the minimum distance to the nearest neighbor
        distances[i] = minimum(dist_to_others)
    end
    
    # Return the mean nearest neighbor distance
    return mean(distances)
end