
# Granular Material Analysis (GranMA)

GranMA is an advanced simulation and data analysis package for the [Julia programming language](https://julialang.org/), designed to provide insights into granular material dynamics. It is compatible with Windows, Linux, and MacOS platforms.

The GranMA ecosystem includes a variety of functionalities:

- **Data Simulation:** GranMA allows for the simulation of granular materials under various conditions, providing detailed insights into their behavior.
- **Data Analysis:** Analyze simulation data to extract meaningful patterns and trends, including frequency, attenuation, and wave number changes.
- **Visualization:** Create detailed plots to visualize the results, aiding in the understanding of complex dynamics.


## Creating Position Vector Structures

The majoritiy of the time we'll spend are with vectors that are either 2D or 3D. We can save significant computational overhead by using `SVector`'s in to define the spatial component of particles.

```julia
using StaticArrays

struct Vec2D{T} <: FieldVector{2,T}
    x::T
    y::T
end
```

We define a specific structure for our position vectors called `Vec2D` that is parameterized be a type `T`.  The type T could be anything, like Float64, Int, etc. The struct will contain two fields, x and y, which will both have the type T. The `<:` symbol indicates that Vec2D is a subtype of `FieldVector{2,T}`, which is an abstract type from `StaticArrays.jl`. We do this so that our `Vec2D` inherits  all the functionalities of a FieldVector, such as element access, iteration, and mathematical operations.

In practice, we simply feed the structure `Vec2D` a type `T` and an `x` and `y` that are type `T` vector is defined.

```julia
julia> Vec2D{Float64}(1,2)
2-element Vec2D{Float64} with indices SOneTo(2):
 1.0
 2.0
```
## Creating Random Vectors

Next, for the sake of making examples, let's create a functions that creates random vectors using the `Vec2D` structure.

```julia
function VecRandom(::Type{VecType},range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    VecRandomOut = VecType(range[begin] + rand(T)*(range[end]-range[begin]) for dim in 1:dimensions)
    return VecRandomOut
end
```


The arguments for our function `VecRandom` are the structure we wan the output vector `VecRandomOut` to be, and the range of the points in that vector. The `where VecType` introduces a type variable `VecType` that can be used throughout the function to refer to the specific type passed as the first argument. This allows us to use different vector types later on, for example 3D.

Here's an example of creating a random vector with this function.

```julia
julia> VecRandom(Vec2D{Float64}, [0,1])
2-element Vec2D{Float64} with indices SOneTo(2):
 0.07920192859937936
 0.8484353399165571
```

## Forces

Now that we have a computationally efficient way to define positions and a way to generate them, we can move on to defining the forces between these points. Let's start with a Hookean spring force law,

$$ F= -k\delta$$

where $\delta$ is the overlap


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
