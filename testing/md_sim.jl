using StaticArrays
using LinearAlgebra: norm

struct Pos2D{T} <: FieldVector{2,T}
    x::T
    y::T
end

struct Pos3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

struct Particle{VecType}
    position::VecType
    diameter::Float64
end

# Function to generate a random position within a given range
function PosRandom(::Type{VecType}, range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    return VecType(range[1] + rand(T)*(range[2]-range[1]) for _ in 1:dimensions)
end

# Function to generate a random particle with random position and diameter
function ParticleRandom(::Type{VecType}, PosRange, DiamRange) where VecType 
    position = PosRandom(VecType, PosRange)
    diameter = DiamRange[1] + rand()*(DiamRange[2] - DiamRange[1])
    return Particle(position, diameter)
end
