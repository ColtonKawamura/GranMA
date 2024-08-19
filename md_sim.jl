using StaticArrays
using LinearAlgebra:norm

struct Vec2D{T} <: FieldVector{3,T}
    x::T
    y::T
    r::T
end

function VecRandom(::Type{VecType},range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    radius = 
    VecRandomOut = VecType(range[begin] + rand(T)*(range[end]-range[begin]) for dim in 1:dimensions)
    return VecRandomOut
end