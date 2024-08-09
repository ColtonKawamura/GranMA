
# Granular Material Analysis (GranMA)

GranMA is an advanced simulation and data analysis package for the [Julia programming language](https://julialang.org/), designed to provide insights into granular material dynamics. It is compatible with Windows, Linux, and MacOS platforms.

The GranMA ecosystem includes a variety of functionalities:

- **Data Simulation:** GranMA allows for the simulation of granular materials under various conditions, providing detailed insights into their behavior.
- **Data Analysis:** Analyze simulation data to extract meaningful patterns and trends, including frequency, attenuation, and wave number changes.
- **Visualization:** Create detailed plots to visualize the results, aiding in the understanding of complex dynamics.


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

```julia
# Clone the repository locally
git clone https://github.com/coltonkawamura/GranMA.git

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

## Julia on HPC
Julia is aleady installed on the HPC, but you'll need to ensure all library depedences for GranMA are installed. Here's how to do it:

After logging into the HPC via ssh, start an interactive session:

```bash
salloc
```

The check the modules availible

```bash
module avail
```

This will bring up all the modules that are availible on the HPC. Look for Julia. At the time on writing this, you might see

```bash
lang/julia/1.9.3   
```

Whatever version is on the HPC is what you'll want to load for the next command you'll put in to load Julia:

```bash
module load lang/julia/1.9.3
```

or whatever version you see.

After this command you can start a julia session with:

```bash
julia
```
Once in the environmnent you can see what packages you already have and load the ones you don't have, for example,

```julia
using MATLAB
```
THis is a good example, because the MATLAB.jl package needs some more work for us to use it compared to other packages. If you followed the steps to install it, you probably got an error. Let's go through the steps to get.

If you're stil in julia, exit out of it to go back to your interactive shell session,

```julia
exit()
```
Then load up the matlab module that appeared in the list of availible modules. Right now for me its

```bash
module load app/matlab/R2023b
```
Next, let's get the path to that module

```bash
which matlab
```

Now add that root to MATLAB_ROOT. For me that path was,

```bash
export MATLAB_ROOT="/share/apps/matlab/R2023b"
```
Now go back into julia and you'll be able to load the MATLAB.jl package!
## Examples

Here are a few examples to get you started with GranMA.

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
