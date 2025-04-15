# Data Analysis

This guide details the data analysis capabilities of GranMA, which are primarily implemented in Julia.

## Overview

The Julia component of GranMA provides powerful tools for:
- Processing raw simulation outputs from MATLAB
- Filtering and selecting data based on physical parameters
- Computing wave propagation metrics (attenuation, dispersion)
- Generating publication-quality visualizations
- Performing advanced statistical analysis

## Core Data Structures

### 2D Simulation Data

The `file_data` structure contains processed data from 2D simulations:

```julia
struct file_data
    pressure::Float64
    omega::Float64
    gamma::Float64
    # Other fields...
    amplitude_vector_x::Vector{Float64}
    amplitude_vector_y::Vector{Float64}
    unwrapped_phase_vector_x::Vector{Float64}
    unwrapped_phase_vector_y::Vector{Float64}
    # Position data...
end
```

### 3D Simulation Data

The `data3d` structure contains processed data from 3D simulations:

```julia
struct data3d
    pressure::Float64
    omega::Float64
    gamma::Float64
    # Other fields...
    amplitude_vector_x::Vector{Float64}
    amplitude_vector_y::Vector{Float64}
    amplitude_vector_z::Vector{Float64}
    unwrapped_phase_vector_x::Vector{Float64}
    unwrapped_phase_vector_y::Vector{Float64}
    unwrapped_phase_vector_z::Vector{Float64}
    # Initial position coordinates...
end
```

## Data Processing Pipeline

### 1. Converting Raw Data

Convert raw MATLAB output files (.mat) to structured Julia data:

```julia
# Process 2D simulation results
simulation_data = crunch("path/to/simulation_output/")

# Process 3D simulation results
simulation_data_3d = crunch3d("path/to/simulation_output/")
```

### 2. Saving & Loading Processed Data

```julia
# Save processed data
save_data(simulation_data, "out/processed/my_data")

# Load data for analysis
data = load_data("out/processed/my_data")
```

### 3. Filtering Data

#### Parameter-Based Filtering

Filter data based on simulation parameters:

```julia
# Filter by pressure, omega, and gamma
filtered_data = FilterData(simulation_data, 
    0.001, :pressure,  # Pressure value, field name
    0.1, :omega,       # Omega value, field name
    0.1, :gamma)       # Gamma value, field name
```

#### Spatial Filtering (3D)

Filter particles based on initial spatial positions:

```julia
# Filter by Y and Z position
filtered_data_3d = FilterData3d(simulation_data_3d,
    0.001, :pressure,
    0.1, :omega,
    0.1, :gamma,
    1, :seed,
    [5, 10], "y",    # Y range [min, max]
    [0, 2], "z")     # Z range [min, max]
```

## Analysis Functions

### Wave Propagation Analysis

```julia
# Calculate wave speed from phase data
wave_speed = calculateWaveSpeed(filtered_data)

# Calculate attenuation from amplitude decay
attenuation = calculateAttenuation(filtered_data)

# Compute dispersion relation
dispersion = computeDispersionRelation(filtered_data)
```

### Statistical Analysis

```julia
# Calculate coordination number statistics
coordination_stats = calculateCoordinationStats(filtered_data)

# Analyze force chain network
force_network = analyzeForceNetwork(filtered_data)

# Compute spatial correlation functions
correlation = spatialCorrelation(filtered_data)
```

## Visualization

### Basic Plots

```julia
# Plot mean field (amplitude and phase vs. distance)
plotMeanField(filtered_data)

# Plot attenuation vs. omega*gamma
plotAttenuationVsOmegaGamma(simulation_data)

# Plot wave speed vs. pressure
plotWaveSpeedVsPressure(simulation_data)
```

### Advanced Visualizations

```julia
# Generate space-time plot of wave propagation
plotSpaceTime(filtered_data)

# Visualize force chains
visualizeForceChains(filtered_data)

# Create animations of wave propagation
createWavePropagationAnimation(filtered_data)
```

### Customizing Plots

```julia
# Example of customizing plot appearance
plotMeanField(filtered_data, 
    title="Wave Propagation in Granular Medium",
    xlabel="Distance from Source",
    ylabel="Amplitude",
    color=:blue,
    legend=true,
    grid=true)
```

## Special Analysis Techniques

### Frequency Domain Analysis

```julia
# Compute and plot frequency spectrum
plotFrequencySpectrum(filtered_data)

# Apply bandpass filtering
filtered_signal = applyBandpassFilter(signal, low_cutoff, high_cutoff)
```

### Scaling Analysis

```julia
# Test power law scaling relationships
exponents = fitPowerLaw(x_data, y_data)

# Plot scaling relationships
plotScalingRelationships(simulation_data)
```

### Comparative Analysis

```julia
# Compare different parameter regimes
compareParameterRegimes(data_set_1, data_set_2, parameter_name)

# Compare to theoretical predictions
compareToTheory(experimental_data, theoretical_model)
```

## Batch Processing

For processing multiple datasets:

```julia
# Process all datasets in a directory
function processAllDatasets(root_dir)
    datasets = []
    for subdir in readdir(root_dir)
        path = joinpath(root_dir, subdir)
        if isdir(path)
            data = crunch(path)
            push!(datasets, data)
        end
    end
    return datasets
end

# Example usage
all_datasets = processAllDatasets("out/simulations/")
```

## Integration with Other Tools

### Exporting to CSV

```julia
# Export data for use in other software
function exportToCsv(data, filename)
    open(filename, "w") do io
        # Write header
        write(io, "pressure,omega,gamma,attenuation,wavespeed\n")
        
        # Write data rows
        for item in data
            row = "$(item.pressure),$(item.omega),$(item.gamma),$(item.attenuation),$(item.wavespeed)\n"
            write(io, row)
        end
    end
end

# Example usage
exportToCsv(filtered_data, "results.csv")
```

### Integration with Python

For workflows that require Python libraries:

```julia
# Using PyCall to interface with Python
using PyCall

# Import Python modules
np = pyimport("numpy")
scipy_signal = pyimport("scipy.signal")

# Example: Use scikit-learn for machine learning
sklearn = pyimport("sklearn.cluster")
kmeans = sklearn.KMeans(n_clusters=3)
# Process data with Python ML tools
```

## Common Analysis Workflows

### Wave Attenuation Analysis

```julia
# Complete workflow for analyzing attenuation vs. omega*gamma
function analyzeAttenuation(data_path)
    # Load data
    data = load_data(data_path)
    
    # Group by pressure
    pressures = unique([d.pressure for d in data])
    
    results = []
    for p in pressures
        # Filter by pressure
        filtered = FilterData(data, p, :pressure)
        
        # Calculate attenuation vs. omega*gamma
        attenuation_data = calculateAttenuationVsOmegaGamma(filtered)
        push!(results, attenuation_data)
        
        # Plot for this pressure
        plotAttenuationVsOmegaGamma(filtered, 
            title="Attenuation vs. omega*gamma (P=$(p))")
    end
    
    return results
end
```

## References

- See the [Code Documentation](./Code-Documentation.md) for detailed function references
- Check Julia's documentation for package-specific functions using `?function_name`