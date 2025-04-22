module EntityTools

using HDF5
using CairoMakie
using Printf
using StatsBase
using Unitful
using Distributions
using StatsBase

include("calc/phase.jl")
include("calc/spectra.jl")
include("debug/parse_timing.jl")

export phase_map, #
        spectrum,
        parse_timing

end
