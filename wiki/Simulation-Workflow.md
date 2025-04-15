# Simulation Workflow

This guide outlines the end-to-end workflow for running granular media acoustic simulations and analyzing the results using the GranMA framework.

## Overview

The GranMA simulation workflow consists of three main stages:

1. **Preparation**: Setting up simulation parameters and generating initial particle configurations
2. **Simulation**: Running the MATLAB-based molecular dynamics simulations
3. **Analysis**: Processing and visualizing results using Julia

## 1. Preparation

### Parameter Configuration

Create a configuration file or script to define your simulation parameters:

```matlab
% Example configuration for a 2D simulation
config = struct();

% Physical parameters
config.width = 40;                     % System width
config.pressure = 0.001;               % Confining pressure
config.spring_constant = 100;          % Spring constant for particle interactions
config.gamma = 0.1;                    % Damping coefficient
config.omega = 0.1;                    % Angular frequency of driving

% Numerical parameters
config.dt = 1e-3;                      % Time step
config.totalTime = 100;                % Total simulation time
config.outputInterval = 0.1;           % Time interval for output
config.seed = 1;                       % Random seed for initial configuration

% Wave parameters
config.amplitude = 0.001;              % Amplitude of driving oscillation
config.driving_type = 'compression';   % 'compression' or 'shear'

% Save configuration
save('in/simulation_config.mat', '-struct', 'config');
```

### Initial Configuration Generation

Generate the initial particle configuration:

```matlab
% Generate initial particle configuration
generateParticleConfiguration(config);
```

## 2. Simulation

### Running a Simulation

To run a single simulation:

```matlab
% Load configuration
config = load('in/simulation_config.mat');

% Run simulation
runGranularSimulation(config);

% Output will be saved in 'out/simulation_2d/' or 'out/simulation_3d/'
```

### Batch Simulations

For parameter sweeps or multiple realizations:

```matlab
% Define parameter ranges
pressures = [0.001, 0.01, 0.1];
omegas = [0.05, 0.1, 0.2, 0.4];
gammas = [0.05, 0.1, 0.2];
seeds = 1:10;

% Run parameter sweep
runParameterSweep(pressures, omegas, gammas, seeds);
```

### HPC Submission

For large simulations on computing clusters:

```bash
# Submit a batch job using SLURM
sbatch hpc_scripts/run_simulation_batch.sh
```

## 3. Analysis

### Data Processing in Julia

After simulations are complete, use Julia for data analysis:

```julia
# Start Julia and activate the project
using Pkg
Pkg.activate("myenv")

# Load the GranMA analysis modules
include("src/GranMA.jl")
using .GranMA

# Process raw simulation data
simulation_data = crunch("out/simulation_2d/")
save_data(simulation_data, "out/processed/data_2d")

# Or for 3D data
simulation_data_3d = crunch3d("out/simulation_3d/")
saveData3d(simulation_data_3d, "out/processed/data_3d")
```

### Data Filtering

Filter processed data for specific analysis:

```julia
# Load processed data
simulation_data = load_data("out/processed/data_2d")

# Filter for specific parameters
filtered_data = FilterData(simulation_data, 
    0.001, :pressure, 
    0.1, :omega, 
    0.1, :gamma)

# For 3D data with spatial filtering
filtered_data_3d = FilterData3d(simulation_data_3d,
    0.001, :pressure, 
    0.1, :omega, 
    0.1, :gamma,
    1, :seed,
    [5, 10], "y",    # Y position range
    [0, 2], "z")     # Z position range
```

### Visualization

Generate plots and visualizations:

```julia
# Generate mean field plots
plotMeanField(filtered_data)

# For 3D data
plotMeanField3d(filtered_data_3d, "y")  # Transverse axis = y

# Generate attenuation vs omega*gamma plots
plotAttenuationVsOmegaGamma(simulation_data)

# Save plots
savefig("figures/meanField_compression_2d.pdf")
```

## Output Analysis

### Key Metrics

The main output metrics include:

- **Wave Speed**: Extracted from phase measurements
- **Attenuation**: Calculated from amplitude decay with distance
- **Dispersion Relations**: Relationship between frequency and wave number
- **Energy Dissipation**: Quantification of energy loss mechanisms

### Visualization Outputs

Common visualization outputs:

1. **Mean Field Plots**: Amplitude and phase vs. distance
2. **Attenuation Plots**: Alpha/omega vs. omega*gamma
3. **Wave Propagation Animations**: Time evolution of wave propagation
4. **Force Chain Visualizations**: Spatial distribution of contact forces

## Advanced Workflows

### Multi-Scale Analysis

For connecting micro- and macro-scale behaviors:

```julia
# Analyze local coordination number effects
analyzeCoordinationNumberEffects(filtered_data)

# Analyze force chain structure
analyzeForceChains(filtered_data)
```

### Machine Learning Integration

For extracting patterns from simulation results:

```julia
# Extract features from simulation data
features = extractSimulationFeatures(simulation_data)

# Train prediction model
trainWavePropagationModel(features, targets)
```

## References

For detailed function documentation:

- See the [Code Documentation](./Code-Documentation.md) wiki page
- Check function-specific help in MATLAB using `help function_name`
- Check function-specific help in Julia using `?function_name`