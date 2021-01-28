# define orbit by its length and minimum value
mutable struct Orbit{T<:AbstractFloat}
    β::T
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

function normalise_basins!(orbits::Array{Orbit,1})
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
    print(io,"Orbit{$T,β=$(repr(o.β))}(length=$(o.length), min=$(repr(o.min)), basin=$(o.b))")
end