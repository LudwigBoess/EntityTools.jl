abstract type DataInfo end

struct FieldInfo <: DataInfo
    fields::Vector{String}
end

struct ParticleInfo <: DataInfo
    species::Vector{Int64}
    properties::Vector{String}
end

struct SpectrumInfo <: DataInfo
    species::Vector{Int64}
end

"""
    struct EntityData 
        path::String
        steps::Vector{Int64}
        times::Vector{Float64}
        settings::Dict{String, Any}
        fields::Union{FieldInfo, Nothing}
        particles::Union{ParticleInfo, Nothing}
        spectra::Union{SpectrumInfo, Nothing}
    end


"""
struct EntityData 
    path::String
    steps::Vector{Int64}
    times::Vector{Float64}
    settings::Dict{String, Any}
    fields::Union{FieldInfo, Nothing}
    particles::Union{ParticleInfo, Nothing}
    spectra::Union{SpectrumInfo, Nothing}
end

function read_times(sim_path, subfolder, read_steps=false)
    # find all available files in the folder
    files = readdir(joinpath(sim_path, subfolder))
    Nfiles = length(files)
    # read written steps
    i_steps = Vector{Int64}(undef, Nfiles)
    for i in 1:Nfiles 
        i_steps[i] = parse(Int64, split(files[i], ".")[2])
    end

    # read times
    times = Vector{Float64}(undef, Nfiles)
    for i = 1:Nfiles
        file = adios_open_serial(joinpath(sim_path, subfolder, "$subfolder.$(@sprintf("%08i", i_steps[i])).bp"), mode_readRandomAccess)
        times[i] = adios_load(file, "Time", 0)
        close(file)
    end

    if read_steps
        return times, i_steps
    else
        return times
    end
end

function read_field_info(sim_path; times=nothing, i_steps=nothing)

    times, i_steps = read_times(sim_path, "fields", true)

    file = adios_open_serial(joinpath(sim_path, "fields", "fields.$(@sprintf("%08i", i_steps[1])).bp"), mode_readRandomAccess)
    variables = adios_all_variable_names(file)
    out_fields = variables[startswith.(variables, "f")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end

    return FieldInfo(out_fields)
end

function read_particle_info(sim_path)

    file = adios_open_serial(joinpath(sim_path, "particles", "particles.$(@sprintf("%08i", 1)).bp"), mode_readRandomAccess)
    variables = adios_all_variable_names(file)
    out_fields = variables[startswith.(variables, "p")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end

    # X1 is always written
    species = [parse(Int, sp[end]) for sp in out_fields[startswith.(out_fields, "X1")]]

    # get unique properties
    properties = unique([pr[1:end-2] for pr in out_fields])

    return ParticleInfo(species, properties)
end

function read_spectrum_info(sim_path)

    file = adios_open_serial(joinpath(sim_path, "spectrum", "spectrum.$(@sprintf("%08i", 1)).bp"), mode_readRandomAccess)
    variables = adios_all_variable_names(file)
    out_fields = variables[startswith.(variables, "s")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end

    # N is always written
    species = [parse(Int, sp[end]) for sp in out_fields[startswith.(out_fields, "N")]]

    return SpectrumInfo(species)
end

function read_run_settings(sim_path, subfolder)
    file = adios_open_serial(joinpath(sim_path, subfolder, "$subfolder.$(@sprintf("%08i", 1)).bp"), mode_readRandomAccess)
    attrs = adios_all_attribute_names(file)
    settings = Dict{String, Any}()
    for attr in attrs
        settings[attr] = adios_load(file, attr, 0)
    end
    close(file)
    return settings
end

"""
    EntityData(sim_path::String)

Reads an overview of the simulation into an `EntityData` struct.
"""
function EntityData(sim_path::String)

    subfolders = readdir(sim_path)
    if iszero(length(subfolders))
        error("No output subfolders found in $sim_path")
    end

    # reading output files and times
    times, i_steps = read_times(sim_path, subfolders[1], true)

    field_info = nothing
    if "fields" in subfolders
        field_info = read_field_info(sim_path)
    end
    particle_info = nothing
    if "particles" in subfolders
        particle_info = read_particle_info(sim_path)
    end
    spectrum_info = nothing
    if "spectrum" in subfolders
        spectrum_info = read_spectrum_info(sim_path)
    end

    # read the settings of the simulation
    runsettings = Dict("dummy" => 0) #read_run_settings(sim_path, subfolders[1])

    return EntityData(sim_path, i_steps, times, runsettings, field_info, particle_info, spectrum_info)
end



