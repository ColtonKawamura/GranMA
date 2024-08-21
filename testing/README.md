
# Granular Material Analysis (GranMA)

## Creating Position Vector Structures

We're only dealing with with vectors that are either 2D or 3D for position. So, we can save significant computational overhead by using `FieldVector`'s in to define the spatial component of particles.

```julia
using StaticArrays

struct Pos2D{T} <: FieldVector{2,T}
    x::T
    y::T
end
```

We define a specific structure for our position vectors called `Pos2D` that is parameterized be a type `T`.  The type T could be anything, like Float64, Int, etc. The struct will contain two fields, x and y, which will both have the type T. The `<:` symbol indicates that Vec2D is a subtype of `FieldVector{2,T}`, which is an abstract type from `StaticArrays.jl`. We do this so that our `Pos2D` inherits  all the functionalities of a FieldVector, such as element access, iteration, and mathematical operations.

In practice, we simply feed the structure `Pos2D` a type `T` and an `x` and `y` that are type `T` vector is defined.

```julia
julia> Pos2D{Float64}(1,2)
2-element Pos2D{Float64} with indices SOneTo(2):
 1.0
 2.0
```

Next, 
## Creating Random Particles

Next, for the sake of making examples, let's create a functions that creates random vectors using the `Vec2D` structure.

```julia
function PosRandom(::Type{VecType}, range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    return VecType(range[1] + rand(T)*(range[2]-range[1]) for _ in 1:dimensions)
end
```


The arguments for our function `PosRandom` are the structure we want the output vector `VecType` to be, and the range of the points in that vector. The `where VecType` introduces a type variable `VecType` that can be used throughout the function to refer to the specific type passed as the first argument. This allows us to use different vector types later on, for example 3D.

Here's an example of creating a random vector with this function.

```julia
julia> PosRandom(Pos2D{Float64}, [0,1])
2-element Pos2D{Float64} with indices SOneTo(2):
 0.28591930991040737
 0.7042845413535875
```

Now that we have a position for the particle, we'll want to keep track of other information of the particle as well, such as the diameter. We can do this by creating a `stuct` that encompases all that informaiton.

```julia
struct Particle{VecType}
    position::VecType
    diameter::Float64
end
```
The structure `Particle` will take on the type `VecType` which is based on the position vector and has fields of `position` and `diameter`. Later, we can include more information if we want each particle do have different masses or spring contants, etc.

Again, to make things easier, let's make a function that will generate particles with random positions and random diamters given a range for each.

```julia
function ParticleRandom(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end
```

This has a similar flow as `PosRandom` function, except it has another argument for the range of the diameter. Let's make a random particle using this function.

```julia
julia> ExampleParticle = ParticleRandom(Pos2D{Float64}, [0,1], [0,1])
Particle{Pos2D{Float64}}([0.7197480175465489, 0.36301736551172714], 0.0846518956414718)
```

We can access the position of the `ExampleParticle` with

```julia
julia> ExampleParticle.position
2-element Pos2D{Float64} with indices SOneTo(2):
 0.7197480175465489
 0.36301736551172714
 ```
 and the diameter,

 ```julia
 julia> ExampleParticle.diameter
0.0846518956414718
 ```
## Forces

Now that we have a computationally efficient way to define and create particles, we can move on to defining the forces between these points. Let's start with a Hookean spring force law where $\delta$ is the overlap. For particles that have no overlap ($\delta = 0$), there is no force. If the distance between particles $i$ and $j$ is given by:

$$
r_{i,j} = \sqrt{(x_j - x_i)^2 + (y_j - y_i)^2}
$$

where $x_i, y_i$ are the positions of particle $i$, and $x_j, y_j$ are the positions of particle $j$, then the potential energy function can be defined as:

$$
U(r) = k\left[(r_i + r_j)^2 - r_{i,j} \right]^2 \quad \text{for} \quad r_{i,j} \leq (r_i + r_j)
$$

and

$$
U(r) = 0.0 \quad \text{for} \quad r_{i,j} > (r_i + r_j)
$$

where $r_i$ and $r_j$ are the radii of particles $i$ and $j$ respectively. Computationally, this is

```julia
function U_Hooke(x::T,y::T,cutoff) where T
    delta = y - x
    d = norm(Δv)
    if d > cutoff
        fₓ = zero(T)
    else
        fₓ = 2*(d - cutoff)*(Δv/d)
    end
    return fₓ
end
```

The force exerted on the particles due to this potential is given by:

$$
\vec{F}(\vec{x}_i, \vec{x}_j, k) = k_i k_j \left[r^2 - (r_i + r_j)^2 \right]
$$

for which the corresponding forces are:

$$
\vec{f}_j = -\vec{f}_i
$$


```julia
function Hooke(x::T,y::T,cutoff) where T
    delta = y - x
    d = norm(Δv)
    if d > cutoff
        fₓ = zero(T)
    else
        fₓ = 2*(d - cutoff)*(Δv/d)
    end
    return fₓ
end
```

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
