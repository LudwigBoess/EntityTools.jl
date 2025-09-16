abstract type DataInfo end

struct FieldInfo <: DataInfo
    steps::Vector{Int64}
    times::Vector{Float64}
    fields::Vector{String}
    extent_min::Vector{Float64}
    extent_max::Vector{Float64}
end

struct ParticleInfo <: DataInfo
    steps::Vector{Int64}
    times::Vector{Float64}
    species::Vector{Int64}
    masses::Vector{Float64}
    properties::Vector{String}
end

struct SpectrumInfo <: DataInfo
    steps::Vector{Int64}
    times::Vector{Float64}
    species::Vector{Int64}
end


struct EntityData 
    path::String
    steps::Vector{Int64}
    times::Vector{Float64}
    settings::Dict{String, Any}
    fields::Union{FieldInfo, Nothing}
    particles::Union{ParticleInfo, Nothing}
    spectra::Union{SpectrumInfo, Nothing}
end