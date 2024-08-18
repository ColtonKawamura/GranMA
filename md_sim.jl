using StaticArrays

struct Vec2D{T} <: FieldVector{2,T}
    x::T
    y::T
end

function VecRandom(::Type{VecType},range) where VecType 
    T = eltype(VecType)
    dimensions = length(VecType)
    VecRandomOut = VecType(range[begin] + rand(T)*(range[end]-range[begin]) for dim in 1:dimensions)
    return VecRandomOut
end