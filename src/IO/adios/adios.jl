abstract type DataInfo end


"""
    struct FieldInfo <: DataInfo
        fields::Vector{String}
    end

Struct to store info on field output.
"""
struct FieldInfo <: DataInfo
    fields::Vector{String}
end

"""
    struct ParticleInfo <: DataInfo
        species::Vector{Int64}
        properties::Vector{String}
    end

Struct to store info on particle output.
"""
struct ParticleInfo <: DataInfo
    species::Vector{Int64}
    properties::Vector{String}
end

"""
    struct SpectrumInfo <: DataInfo
        species::Vector{Int64}
    end

Struct to store info on output spectra.
"""
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

Struct to store metadata of an Entity simulation.
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
    read_times(sim_path, subfolder, read_steps=false)

Finds all available output times for a given `subfolder`.
If `read_steps=true` is set it also returns an array of the output steps constructed from file names.
"""
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
    find_closest_time(run_info::EntityData, t::Real)

Finds the time in the available outputs that is closest to the requested one
"""
function find_closest_time(run_info::EntityData, t::Real)
    return run_info.times[argmin(abs.(run_info.times .- t))]
end

"""
    find_idxs(run_info, i_step, t, i, reduction_function, verbose)

Find requested output file(s) based on single/array of steps `i_step`, time `t` or file number `i`.
If a `reduction_function` is provided, options are interpreted as arrays.
"""
function find_idxs(run_info::EntityData, i_step, t, i, reduction_function, verbose)

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
            error("You can only read one file in this mode. To get the properties over multiple files please provide a `reduction_function``.")
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

"""
    select_first_file(sim_path, subfolder)

Selects the first available file in a `subfolder` or throws an error 
"""
function select_first_file(sim_path, subfolder)

    # check if there are files present
    folder_files = readdir(joinpath(sim_path, subfolder))
    if isempty(folder_files)
        error("No files present in $(subfolder)!")
    end
    # open first available file
    file = adios_open_serial(joinpath(sim_path, subfolder, folder_files[1]), mode_readRandomAccess)

    return file
end

"""
    read_field_info(sim_path)

Reads the info on what fields are stored in the `particles` files.
"""
function read_field_info(sim_path)

    # open first available file
    file = select_first_file(sim_path, "fields")
    # read all available output fields
    variables = adios_all_variable_names(file)
    # field fields start with `f`
    out_fields = variables[startswith.(variables, "f")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end

    return FieldInfo(out_fields)
end

"""
    read_particle_info(sim_path)

Reads the info on what fields are stored in the `particles` files.
"""
function read_particle_info(sim_path)

    # open first available file
    file = select_first_file(sim_path, "particles")
    # read all available output fields
    variables = adios_all_variable_names(file)
    # particle fields start with `p`
    out_fields = variables[startswith.(variables, "p")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end
    # X1 is written for every species in all dimensions
    species = [parse(Int, sp[end]) for sp in out_fields[startswith.(out_fields, "X1")]]
    # get unique properties
    properties = unique([pr[1:end-2] for pr in out_fields])

    return ParticleInfo(species, properties)
end


"""
    read_spectrum_info(sim_path)

Reads the info on what fields are stored in the `spectrum` files.
"""
function read_spectrum_info(sim_path)

    # open first available file
    file = select_first_file(sim_path, "spectra")
    # read all available output fields
    variables = adios_all_variable_names(file)
    # spectrum fields start with `s`
    out_fields = variables[startswith.(variables, "s")]
    for i = 1:length(out_fields)
        out_fields[i] = out_fields[i][2:end]
    end

    # N is always written
    species = [parse(Int, sp[end]) for sp in out_fields[startswith.(out_fields, "N")]]

    return SpectrumInfo(species)
end

"""
    read_run_settings(sim_path, subfolder)

Read the general settings of a run from the first file in a subfolder.
"""
function read_run_settings(sim_path, subfolder)
    
    # open first available file
    file = select_first_file(sim_path, subfolder)
    # read run setting names
    attrs = adios_all_attribute_names(file)
    # allocate storage dict
    settings = Dict{String, Any}()
    # read all settings into the dict
    for attr in attrs
        settings[attr] = adios_attribute_data(file, attr)
    end
    close(file)
    return settings
end

"""
    EntityData(sim_path::String)

Reads an overview of the simulation into an `EntityData` struct.
"""
function EntityData(sim_path::String)

    if !isdir(sim_path)
        error("$sim_path not present!")
    end

    # check which output folders are present
    all_files = readdir(sim_path)

    # filter out only the directories
    subfolders = filter(entry -> isdir(joinpath(sim_path, entry)), all_files)

    # make sure to not try to read a plots folder
    filter!(x -> x != "Plots", subfolders)
    filter!(x -> x != "plots", subfolders)

    # check if any folders are present
    if iszero(length(subfolders))
        error("No output subfolders found in $sim_path")
    end

    # reading output files and times
    times, i_steps = read_times(sim_path, subfolders[1], true)

    # if fields are written, read their info
    field_info = nothing
    if "fields" in subfolders
        field_info = read_field_info(sim_path)
    end
    
    # if particles are written, read their info
    particle_info = nothing
    if "particles" in subfolders
        particle_info = read_particle_info(sim_path)
    end

    # if spectra are written, read their info
    spectrum_info = nothing
    if "spectra" in subfolders
        spectrum_info = read_spectrum_info(sim_path)
    end

    # read the settings of the simulation
    runsettings = read_run_settings(sim_path, subfolders[1])

    # construct an EntityData struct
    return EntityData(sim_path, i_steps, times, runsettings, field_info, particle_info, spectrum_info)
end