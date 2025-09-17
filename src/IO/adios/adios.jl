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

"""
    find_idxs(run_info, i_step, t, i, reduction_function, verbose)


"""
function find_idxs(run_info, i_step, t, i, reduction_function, verbose)

    if isnothing(i_step) && isnothing(t) && isnothing(i)
        error("Must provide either i_step, time t, or index i to read a field.")
    end
        
    # case for multiple files
    if !isnothing(reduction_function)
        if !isnothing(i_step)
            if length(i_step) == 2
                i = findfirst(isequal(i_step[1]), run_info.steps):findfirst(isequal(i_step[2]), run_info.steps)
            else
                i = [findfirst(isequal(i_step[j]), run_info.steps) for j = 1:length(i_step)]
            end
        end
        if !isnothing(t)
            if length(t) == 2
                i = findall(x -> (x >= t[1]) && (x <= t[2]), run_info.times)
            else
                i = [argmin(abs.(run_info.times .- t[j])) for j = 1:length(t)]
            end
        end

        if i[1] < 1 || i[end] > length(run_info.steps)
            error("Index i=$i is out of bounds. Must be between 1 and $(length(run_info.steps)).")
        end
    else # no reduction function -> single step
        if (!isnothing(i_step) && length(i_step) > 1) || 
            (!isnothing(t) &&  length(t) > 1) || 
            (!isnothing(i) && length(i) > 1)
            error("You can only read one file in this mode. To get the properties over multiple files please provide a reduction_function.")
        end

        if !isnothing(i_step)
            i = findfirst(isequal(i_step), run_info.steps)
            if isnothing(i)
                error("Step $i_step not found in field data.")
            end
        end
        if !isnothing(t)
            i = argmin(abs.(run_info.times .- t))
            if verbose
                @info "Requested time $t: reading closest time $(run_info.times[i]) at step $(run_info.steps[i])."
            end
        end

        if i < 1 || i > length(run_info.steps)
            error("Index i=$i is out of bounds. Must be between 1 and $(length(run_info.steps)).")
        end
    end

    return i

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



