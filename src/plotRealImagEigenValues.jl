using Plots
# Use PlotlyJS backend for interactive, zoomable plots
plotlyjs()
using ColorSchemes
using LinearAlgebra
using LaTeXStrings

"""
    plotRealImagEigenValues(plotData, pressure_list, damping_list; loglog=false, semilogx=false, semilogy=false, scaling=false)

Plot real and imaginary parts of eigenvalues for granular systems.

# Arguments
- `plotData`: Dictionary containing eigenvalues and other properties
- `pressure_list`: Array of pressures to plot
- `damping_list`: Array of damping constants to plot
- `loglog`: Whether to use logarithmic scale for both axes
- `semilogx`: Whether to use logarithmic scale for x-axis
- `semilogy`: Whether to use logarithmic scale for y-axis
- `scaling`: Whether to scale the eigenvalues by pressure and damping

# Example
```julia
plotRealImagEigenValues(outData, [0.1, 0.2], [0.5, 1.0], loglog=true)
```
"""
function plotRealImagEigenValues(plotData, pressure_list, damping_list; 
                               loglog=false, semilogx=false, semilogy=false, scaling=false)
    
    # Create the plot
    p = plot(layout=(1,1), size=(800, 600), margin=10Plots.mm, legend=:outertopright)
    
    # Get appropriate color scheme
    colors = cgrad(ColorSchemes.viridis, length(pressure_list) > 1 ? length(pressure_list) : length(damping_list))
    
    for (i, pressure) in enumerate(pressure_list)
        for (j, damping_constant) in enumerate(damping_list)
            # Find matching data for the given pressure and damping
            # Handle MATLAB imported arrays which might have multiple dimensions
            pressure_array = plotData["pressure"]
            damping_array = plotData["damping"]
            
            # Convert to vectors if they aren't already
            if ndims(pressure_array) > 1
                pressure_vec = vec(pressure_array)
            else
                pressure_vec = pressure_array
            end
            
            if ndims(damping_array) > 1
                damping_vec = vec(damping_array)
            else
                damping_vec = damping_array
            end
            
            # Use approximate matching with tolerance
            tol = 1e-5
            matching_indices = findall(x -> (abs(pressure_vec[x] - pressure) < tol) && 
                                           (abs(damping_vec[x] - damping_constant) < tol),
                                      1:length(pressure_vec))
            
            if isempty(matching_indices)
                # Try with a looser tolerance
                tol = 1e-3
                matching_indices = findall(x -> (abs(pressure_vec[x] - pressure) < tol) && 
                                               (abs(damping_vec[x] - damping_constant) < tol),
                                          1:length(pressure_vec))
                
                if isempty(matching_indices)
                    @warn "No data found for pressure=$pressure and damping=$damping_constant"
                    println("Available pressure values: ", unique(round.(pressure_vec, digits=4)))
                    println("Available damping values: ", unique(round.(damping_vec, digits=4)))
                    continue
                end
            end
            
            idx = matching_indices[1]  # Use the first matching set of data
            
            # Extract eigenvalues - handle different possible formats from MATLAB
            eigenValues = plotData["eigenValues"][idx]
            
            # Check if we need to convert from a matrix to a vector
            if ndims(eigenValues) > 1
                eigenValues = vec(eigenValues)
            end
            
            # Handle complex values
            realEigenValues = real.(eigenValues)
            imagEigenValues = imag.(eigenValues)
            
            # Keep only positive imaginary parts (like in the MATLAB code)
            keepIdx = imagEigenValues .>= 0
            realEigenValues = realEigenValues[keepIdx]
            imagEigenValues = imagEigenValues[keepIdx]
            
            # Get system dimensions
            Lx = plotData["Lx"][idx]
            Ly = plotData["Ly"][idx]
            
            if ndims(Lx) > 0
                Lx = Lx[1]
            end
            if ndims(Ly) > 0
                Ly = Ly[1]
            end
            
            # Determine marker color
            if length(pressure_list) == 1
                marker_color = :blue
            else
                # Scale color based on pressure value
                color_idx = findfirst(x -> x ≈ pressure, sort(unique(pressure_list)))
                marker_color = colors[color_idx]
            end
            
            # Determine marker size
            if length(damping_list) == 1
                markerSize = 4
            else
                # Scale size based on damping
                markerSize = exp(damping_constant / maximum(damping_list)) * 3
            end
            
            # Create label for the plot
            pressureLabel = "\$$(round(pressure, digits=4)), $(round(damping_constant, digits=4))\$"
            pressureValue = pressure
            
            # Plot the data with appropriate scaling
            if scaling
                scatter!(p, 
                    imagEigenValues ./ sqrt(pressureValue), 
                    -realEigenValues ./ (sqrt(pressureValue) .* damping_constant),
                    markersize=markerSize,
                    markerstrokewidth=0,
                    color=marker_color,
                    label=pressureLabel,
                    alpha=0.7
                )
                ylabel!(p, L"Re(\lambda)/(\sqrt{\hat{P}\hat{\gamma}})")
                xlabel!(p, L"Im(\lambda)/\sqrt{\hat{P}}")
            else
                scatter!(p, 
                    imagEigenValues, 
                    -realEigenValues ./ damping_constant,
                    markersize=markerSize,
                    markerstrokewidth=0,
                    color=marker_color,
                    label=pressureLabel,
                    alpha=0.7
                )
                ylabel!(p, L"Re(\lambda)/\hat{\gamma}")
                xlabel!(p, L"Im(\lambda)")
            end
        end
    end
    
    # Apply appropriate axis scaling
    if loglog
        plot!(p, xscale=:log10, yscale=:log10)
    elseif semilogx
        plot!(p, xscale=:log10)
    elseif semilogy
        plot!(p, yscale=:log10)
    end
    
    # Add title with system dimensions
    Lx_val = plotData["Lx"][1]
    Ly_val = plotData["Ly"][1]
    title!(p, "Lx by Ly: $(round(Lx_val, digits=2)) by $(round(Ly_val, digits=2))")
    
    # Set legend title
    plot!(p, legendtitle=L"\hat{P}, \hat{\gamma}")
    
    # Set plot appearance
    plot!(p, grid=true, fontfamily="Computer Modern", guidefontsize=12, legendfontsize=10)
    
    return p
end

"""
    slopeLine(plotType, slope, xrange, offset; TextLocation=nothing)

Add a line with specified slope to an existing plot.

# Arguments
- `plotType`: Type of plot, e.g., "loglog", "linear"
- `slope`: Slope of the line
- `xrange`: Range of x values for the line [xmin, xmax]
- `offset`: Vertical offset of the line
- `TextLocation`: Optional position for slope text label

# Example
```julia
slopeLine("loglog", 1.0, [0.1, 1.0], 0.5, TextLocation=[0.5, 0.5])
```
"""
function slopeLine(plotType, slope, xrange, offset; TextLocation=nothing)
    x_vals = range(xrange[1], xrange[2], length=100)
    
    if plotType == "loglog" || plotType == "semilogx"
        # For log-scale x, use logarithmic spacing
        x_vals = exp.(range(log(xrange[1]), log(xrange[2]), length=100))
    end
    
    if plotType == "loglog" || plotType == "linear"
        y_vals = offset .* (x_vals .^ slope)
    elseif plotType == "semilogx"
        y_vals = offset .+ slope .* log.(x_vals)
    else
        error("Unsupported plot type: $plotType")
    end
    
    # Plot the line
    plot!(x_vals, y_vals, linestyle=:dash, linewidth=2, label=nothing, color=:black)
    
    # Add text label if specified
    if !isnothing(TextLocation)
        tx, ty = TextLocation
        annotate!(tx, ty, text("Slope: $slope", 10, :black))
    end
end

"""
Example usage:

```julia
using MAT

# Load the data
vars = matread("out/2d_damped_eigenStuff/combined/2D_damped_eigenstuff_N1483_40by56_K100_M1.mat")
outData = vars["outData"]

# Plot eigenvalues for specific pressure and damping values 
p = plotRealImagEigenValues(outData, [0.001, 0.005, 0.01], [0.05], loglog=true)
slopeLine("loglog", 1, [0.5, 5], 2)
savefig(p, "eigenvalue_plot.png")
```
```
