module EntityTools

    using CairoMakie
    using Printf
    using StatsBase
    using Unitful
    using Distributions
    using StatsBase
    using ProgressMeter
    using ADIOS2
    using Printf
    using TOML
    using FFTW
    using WriteVTK 
    using Interpolations
    using LinearAlgebra
    using Makie
    using GeometryBasics

    include("IO/adios/adios.jl")
    include("IO/adios/fields.jl")
    include("IO/adios/particles.jl")
    include("IO/adios/spectra.jl")
    include("IO/vtk/write_vtk.jl")
    include("calc/phase.jl")
    include("calc/spectra.jl")
    include("calc/powerspectrum.jl")
    include("debug/parse_timing.jl")
    include("units/structs.jl")
    include("units/setup.jl")
    include("units/length.jl")
    include("units/speed.jl")
    include("units/time.jl")
    include("units/Bfield.jl")
    include("units/temperature.jl")
    include("plotting/polar.jl")

    export phase_map,
            spectrum,
            parse_timing,
            EntityData,
            EntityUnits,
            find_closest_time,
            read_particles,
            read_field,
            read_spectrum,
            bp_to_vtk,
            power_spectrum,
            get_Temp,
            polar_fieldlines!,
            polar_pcolor!,
            polar_contour!

end
