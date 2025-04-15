# Installation Guide

This guide will help you set up your environment to run GranMA simulations and analysis.

## System Requirements

- **Operating Systems**: MacOS, Linux, or Windows
- **RAM**: 8GB minimum, 16GB+ recommended for large simulations
- **Storage**: At least 10GB of free disk space for simulation outputs
- **CPU**: Multi-core processor recommended

## Software Requirements

### MATLAB Setup (for Simulations)

1. **MATLAB Installation**:
   - MATLAB R2019b or later
   - Required Toolboxes:
     - Parallel Computing Toolbox
     - Statistics and Machine Learning Toolbox
     - Signal Processing Toolbox

2. **MATLAB Configuration**:
   ```matlab
   % Add the GranMA simulation paths to your MATLAB environment
   addpath('/path/to/GranMA2025/matlab_scripts');
   savepath;
   ```

### Julia Setup (for Analysis)

1. **Install Julia**:
   - Download and install Julia 1.7+ from [https://julialang.org/downloads/](https://julialang.org/downloads/)
   - Add Julia to your system PATH

2. **Required Julia Packages**:
   ```julia
   # Run this in the Julia REPL
   using Pkg
   Pkg.add([
       "MAT",
       "Plots",
       "StatsBase",
       "LinearAlgebra",
       "FFTW",
       "JLD2",
       "DSP",
       "Interpolations",
       "ProgressMeter"
   ])
   ```

3. **Julia Environment Setup**:
   ```bash
   # Clone the GranMA repository
   git clone https://github.com/yourusername/GranMA2025.git
   cd GranMA2025
   
   # Activate the project environment
   julia -e 'using Pkg; Pkg.activate("myenv"); Pkg.instantiate()'
   ```

## Optional Components

### HPC Environment Setup

For running large-scale simulations on high-performance computing clusters:

1. **SLURM Scripts**: 
   - Sample scripts are provided in the `hpc_scripts/` directory
   - Modify the resource allocation parameters according to your cluster's requirements

2. **Module Loading**:
   ```bash
   # Example for loading required modules on HPC
   module load matlab/R2021a
   module load julia/1.7.2
   ```

### Visualization Tools

For enhanced visualization options:

1. **ParaView**: For 3D visualization of particle configurations
   - Download from [https://www.paraview.org/download/](https://www.paraview.org/download/)

2. **PyVista**: For Python-based visualization workflows
   ```bash
   pip install pyvista matplotlib numpy
   ```

## Verifying Installation

To verify your installation is working correctly:

```bash
# Navigate to the project directory
cd /path/to/GranMA2025

# Run a test simulation in MATLAB
matlab -nodisplay -nosplash -nodesktop -r "run('tests/test_simulation.m'); exit;"

# Run a test analysis in Julia
julia --project=myenv -e "include(\"tests/test_analysis.jl\")"
```

## Troubleshooting

1. **MATLAB Path Issues**:
   - Ensure all required directories are in MATLAB's path
   - Check for path conflicts with other toolboxes

2. **Julia Package Issues**:
   - Try `Pkg.update()` to ensure all packages are on compatible versions
   - For version conflicts, use `Pkg.resolve()` to find a working combination

3. **Performance Issues**:
   - Check that Julia is using multi-threading by setting environment variable: 
     - `export JULIA_NUM_THREADS=4` (adjust to your CPU core count)

For additional help, please refer to the [Troubleshooting](./Troubleshooting.md) page or submit an issue on the GitHub repository.