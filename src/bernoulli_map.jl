import Distributed: @distributed

"""Returns an array of subsequent values from the Bernoulli map with parameter β, starting
from x over N time steps."""
function bernoulli_map(x::T,β::T,N::Int) where T
    oone = one(T)                   # 1 in format T
    xout = Array{T,1}(undef,N+1)    # preallocate
    xout[1] = x                     # store initial condition
    @inbounds for i in 2:N+1
        x = β*x                     # actual Bernoulli map
        x = x >= oone ? x-oone : x  # x mod 1
        xout[i] = x                 # store iteration
    end
    return xout
end

"""Find the Bernoulli orbit length that is reached starting from x with parameter β.
Recursive function that starts testing for orbits over period length N, but enlarges
that testing period if no orbit was found."""     
function orbit_length(  x::T,            # initial condition
                        β::T,            # Bernoulli parameter
                        N::Int) where T  # max period length
    oone = one(T)                          # 1 in format T
    x0 = x                          # new initial condition
    n = 0                           # orbit length (0 = not found yet)
    j = 0                           # iteration counter

    while n == 0 && j < N           # stop when orbit is found (n >0) or when max period is reached
        j += 1
        x = β*x                     # actual Bernoulli map
        x = x >= oone ? x-oone : x  # x mod 1
        n = x0 == x ? j : 0         # check for periodicity
    end

    if n == 0         # in case no orbit was found in N steps, repeat with 5N steps
        n,x = orbit_length(x,β,5N)
    else
        return n,x    # return orbit length and last x that's on the orbit
    end
end

"""Returns the minimum x on the orbit of length N, starting from x (which is assumed to be on the orbit."""
function orbit_minimum( x::T,           # initial condition on the orbit
                        β::T,           # Bernoulli parameter
                        N::Int) where T # the period length of the orbit
    oone = one(T)                   # 1 in format T
    xmin = x                        # store the orbit's minimum in xmin
    for i in 2:N+1
        x = β*x                     # actual Bernoulli map
        x = x >= oone ? x-oone : x  # x mod 1
        xmin = x < xmin ? x : xmin  # always store the smallest x on the orbit
    end
    return xmin
end

"""Returns both the length and minimum of the orbit reached from initial condition x."""
function orbit_length_minimum(  x0::T,              # initial condition on the orbit
                                β::T,               # Bernoulli parameter
                                N::Int=10) where T  # initial period length to test for
    o,x = orbit_length(x0,β,N)    
    x = orbit_minimum(x,β,o)
    return Orbit(β,o,x,Float64(eps(x0)))
end

"""For a given number format T, parameter β and n initial conditions, find orbits of the
Bernoulli map and test for their uniqueness."""
function find_orbits(   ::Type{T},                  # Number format
                        β::Real,                    # Bernoulli parameter
                        n::Int=10000) where T       # n initial conditions

    # pre-allocate empty array of orbits
    orbits = Orbit[]
    uint_type = eval(Symbol("UInt"*repr(sizeof(T)*8)))
    Tβ = T(β)

    tic = time()
    orbits = @distributed (reduce_orbits) for i in 1:n            # for n ICs calculate orbit lengths & x
        orbit_length_minimum(reinterpret(T,uint_type(i)),Tβ,10)    
    end
    
    sort!(orbits)   # from shortest to longest orbit
    normalise_basins!(orbits)
    
    toc = time()
    time_elapsed = readable_secs(toc-tic)
    println("Found $(length(orbits)) orbits in $time_elapsed.")
    return orbits
end

"""For a given number format T, parameter β and n initial conditions, find orbits of the
Bernoulli map and test for their uniqueness."""
function find_orbits_rand(  ::Type{T},                  # Number format
                            β::Real,                    # Bernoulli parameter
                            n::Int=10000) where T       # n initial conditions

    # pre-allocate empty array of orbits
    orbits = Orbit[]
    uint_type = eval(Symbol("UInt"*repr(sizeof(T)*8)))
    Tβ = T(β)

    tic = time()
    orbits = @distributed (reduce_orbits) for i in 1:n            # for n ICs calculate orbit lengths & x
        orbit_length_minimum(T(rand()),Tβ,10)    
    end
    
    sort!(orbits)   # from shortest to longest orbit
    normalise_basins!(orbits)
    
    toc = time()
    time_elapsed = readable_secs(toc-tic)
    println("Found $(length(orbits)) orbits in $time_elapsed.")
    return orbits
end