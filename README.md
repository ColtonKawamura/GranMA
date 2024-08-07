
# Granula Material Analysis (GranMA)

GranMA is an advanced simulation and data analysis package for the [Julia programming language](https://julialang.org/), designed to provide insights into granular material dynamics. It is compatible with Windows, Linux, and MacOS platforms.

The GranMA ecosystem includes a variety of functionalities:

- **Data Simulation:** GranMA allows for the simulation of granular materials under various conditions, providing detailed insights into their behavior.
- **Data Analysis:** Analyze simulation data to extract meaningful patterns and trends, including frequency, attenuation, and wave number changes.
- **Visualization:** Create detailed plots to visualize the results, aiding in the understanding of complex dynamics.


To learn more, we invite you to visit the documentation at [docs.granma.org](http://docs.granma.org/stable/).

## Installation

To get started with GranMA, install the package using Julia's package manager. You can install the latest release of GranMA as follows:

```julia
julia>]
pkg> add GranMA
```

Check the installed version:

```julia
]st GranMA
```

Start using the package:

```julia
using GranMA
```

## Developing GranMA

<details>
  <summary><span style="color:red"> 🔥 Click for more 🔥</span></summary>

GranMA is developed as a modular system, allowing easy customization and extension. To contribute to its development, follow these steps:

```julia
# Clone the repository locally
git clone https://github.com/yourusername/GranMA.git

# Enter the package directory
cd GranMA

# Start the Julia REPL and activate the environment
julia --project=.

# Install dependencies
julia>]
pkg> instantiate

# Run tests
pkg> test
```

Feel free to explore the code, submit issues, and contribute improvements via pull requests!
</details>

## Examples

Here are a few examples to get you started with GranMA. For more detailed information, please refer to [the documentation](http://docs.granma.org/stable/).

### Basic Simulation

```julia
using GranMA

# Set up simulation parameters
parameters = SimulationParameters(
    pressure = 0.01,
    omega = 1.0,
    gamma = 0.1
)

# Run the simulation
result = run_simulation(parameters)

# Visualize results
plot_results(result)
```

### Advanced Analysis

<details>
  <summary>Show Code</summary>

```julia
using GranMA

# Load simulation data
data = load_data("path/to/simulation_output.mat")

# Perform frequency analysis
frequency_results = analyze_frequencies(data)

# Plot frequency spectrum
plot_frequency_spectrum(frequency_results)
```
</details>

## Sponsors

We are grateful to our sponsors for supporting the development of GranMA. Their generous contributions help us continue our mission to provide high-quality tools for granular material analysis.
