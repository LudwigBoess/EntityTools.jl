module EntityTools

    using CairoMakie
    using Printf
    using StatsBase
    using Statistics
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
    using JSON3

    include("IO/adios/adios.jl")
    include("IO/adios/fields.jl")
    include("IO/adios/particles.jl")
    include("IO/adios/spectra.jl")
    include("IO/adios/performance.jl")
    include("IO/vtk/write_vtk.jl")
    include("IO/parse_output.jl")
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
            parse_output,
            SimOutput,
            StepData,
            substep_series,
            total_series,
            active_particles,
            timestep_durations,
            elapsed_times,
            remaining_times,
            sim_times,
            step_numbers,
            species_series,
            summary_table,
            parse_compound_time,
            adios2_throughput,
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
