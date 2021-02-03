# define orbit by its length and minimum value
mutable struct Orbit{T<:AbstractFloat}
    β::Real
    length::Int
    min::T
    b::Float64
end

# two orbits are equivalent when their βs, length and min agree
function Base.:(==)(o1::Orbit{T},o2::Orbit{T}) where T
    (o1.β == o2.β) && (o1.length == o2.length) && (o1.min == o2.min)
end

# two equivalent orbits added adds their basins
function Base.:(+)(o1::Orbit{T},o2::Orbit{T}) where T
    @assert o1 == o2 "$o1 != $o2."
    o1.b += o2.b
    return o1
end

# addition without uniqueness check
function unsafe_add(o1::Orbit{T},o2::Orbit{T}) where T
    o1.b += o2.b
    return o1
end

# return index of first element in x equal to y
function findin(x::Array{T,1},y::T) where T
    @inbounds for i in 1:length(x)
        x[i] == y && return i
    end
end

# orbit array reduction methods for distributed calculation
# check for uniqueness and extend orbit array if so
reduce_orbits(x::Orbit{T},y::Orbit{T}) where T = x == y ? [unsafe_add(x,y)] : [x,y]

function reduce_orbits(x::Orbit{T},y::Array{Orbit{T},1}) where T
    if x in y
        y[findin(y,x)] += x
        return y
    else
        return push!(y,x)
    end
end

function reduce_orbits(x::Array{Orbit{T},1},y::Orbit{T}) where T
    if y in x
        x[findin(x,y)] += y
        return x
    else
        return push!(x,y)
    end
end

function reduce_orbits(x::Array{Orbit{T},1},y::Array{Orbit{T},1}) where T
    for o in x
        y = reduce_orbits(o,y)
    end
    return y
end

function normalise_basins!(orbits::Array{Orbit{T},1}) where T
    s = 0.0
    for o in orbits
        s += o.b
    end
    for o in orbits
        o.b /= s
    end
end
     
# define < for sorting 
function Base.isless(o1::Orbit{T},o2::Orbit{T}) where T
    o1.length < o2.length
end

function Base.show(io::IO, o::Orbit{T}) where T

    l = @sprintf("%6d",o.length)
    m = @sprintf("%16s",repr(o.min))
    b = @sprintf("%22s",repr(o.b))

    print(io,"Orbit{$T,β=$(repr(o.β))}(length=$l, min=$m, basin=$b)")
end