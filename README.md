# Granular Media Acoustics (GranMA)

## Project Overview
A comprehensive computational framework for investigating acoustic wave propagation in granular materials using advanced molecular dynamics (MD) simulations. This research explores how particle-level characteristics and packing arrangements influence macro-scale acoustic properties, with implications for materials science, seismology, and advanced sensing applications.

## Workflow Architecture

The research pipeline consists of three integrated stages:

1. **Packing Generation**: Creating physically realistic granular assemblies with controlled properties
2. **Wave Propagation Simulation**: MATLAB-based molecular dynamics simulations of acoustic wave propagation
3. **Data Analysis & Visualization**: Julia-based post-processing for extracting physical insights

## Technical Implementation

### Simulation Engine (MATLAB)
The core simulation engine leverages MATLAB's robust numerical solvers to efficiently compute:
- Time-evolving particle dynamics with complex contact mechanics
- Hertzian and non-Hertzian contact models with customizable parameters
- Boundary conditions mimicking experimental setups
- Parallelized computation for large particle systems

### Analysis Framework (Julia)
Post-processing is implemented in Julia to capitalize on several key advantages:

- **Performance**: Julia's just-in-time compilation offers near-C speeds while maintaining high-level syntax, enabling rapid processing of large simulation datasets (≥10GB)
- **Scientific Ecosystem**: Julia's purpose-built scientific computing ecosystem provides optimized statistics, signal processing, and visualization libraries
- **Code Clarity**: Julia's expressive syntax allows for concise yet readable implementations of complex analysis algorithms, enhancing reproducibility
- **Future-proofing**: As a modern language explicitly designed for scientific computing, Julia supports the latest computational methods and hardware acceleration techniques

The transition from MATLAB to Julia for analysis represents a strategic technical decision that balances development speed, computational performance, and long-term maintainability for a research project of this scale and complexity.

## Documentation

Comprehensive documentation including theoretical background, usage examples, and API references is available in the [project wiki](https://github.com/ColtonKawamura/GranMA/wiki).

## Research Applications

Current research using this framework focuses on:
- Characterizing wave attenuation in pressure-dependent granular systems
- Investigating the relationship between packing geometry and acoustic anisotropy
- Developing predictive models for wave propagation in complex granular media
