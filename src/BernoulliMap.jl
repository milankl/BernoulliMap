module BernoulliMap

    export Orbit, find_orbits, find_orbits_rand,
        bernoulli_map, orbit_length_minimum, orbit_minimum, orbit_length

    include("utils.jl")
    include("orbits.jl")
    include("bernoulli_map.jl")
end
