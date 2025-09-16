function read_particles_reduced(run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose=true)

    if isnothing(i_step) && isnothing(t) && isnothing(i)
        error("Must provide either i_step, time t, or index i to read a field.")
    end

    if !isnothing(i_step)
        if length(i_step) == 2
            i = findfirst(isequal(i_step[1]), run_info.field_info.i_steps):findfirst(isequal(i_step[2]), run_info.field_info.i_steps)
        else
            i = [findfirst(isequal(i_step[j]), run_info.field_info.i_steps) for j = 1:length(i_step)]
        end
    end
    if !isnothing(t)
        if length(t) == 2
            i = findall(x -> (x >= t[1]) && (x <= t[2]), run_info.field_info.times)
        else
            i = [argmin(abs.(run_info.field_info.times .- t[j])) for j = 1:length(t)]
        end
    end

    if i < 1 || i > length(run_info.field_info.i_steps)
        error("Index i=$i is out of bounds. Must be between 1 and $(length(run_info.field_info.i_steps)).")
    end

    data = Dict(String, Vector{Float64})()
    data["t"] = [run_info.field_info.times[i]]
    data["step"] = [run_info.field_info.steps[i]]

    for j in i
        # construct the file name
        fi = joinpath(RunInfo.path, "particles", "particles.$(@sprintf("%08i", run_info.field_info.i_steps[j])).bp")
        # read the field
        file = adios_open_serial(fi, mode_readRandomAccess)
        for prop in properties
            field_name = "$(prop)_$(species)"
            if !(field_name in run_info.particle_info.properties)
                error("Property $prop for species $species not found in particle data.")
            end
            data[prop][j] = reduction_function(adios_load(file, "p$field_name", 0))
        end    # remember to close the file
        close(file)

    end

    return field
end


function read_particles(run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose=true)

    if isnothing(i_step) && isnothing(t) && isnothing(i)
        error("Must provide either i_step, time t, or index i to read a field.")
    end

    if !isnothing(reduction)
        return read_particles_reduced(run_info, species, properties; 
                    i_step=i_step, t=t, i=i,
                    reduction_function=reduction_function,
                    verbose=verbose)
    end

    if !isnothing(i_step)
        i = findfirst(isequal(i_step), run_info.field_info.i_steps)
        if isnothing(i)
            error("Step $i_step not found in field data.")
        end
    end
    if !isnothing(t)
        i = argmin(abs.(run_info.field_info.times .- t))
        if verbose
            @info "Requested time $t: reading closest time $(run_info.field_info.times[i]) at step $(run_info.field_info.steps[i])."
        end
    end

    if i < 1 || i > length(run_info.field_info.i_steps)
        error("Index i=$i is out of bounds. Must be between 1 and $(length(run_info.field_info.i_steps)).")
    end

    data = Dict(String, Vector{Float64})()
    data["t"] = [run_info.field_info.times[i]]
    data["step"] = [run_info.field_info.steps[i]]

    # construct the file name
    fi = joinpath(RunInfo.path, "particles", "particles.$(@sprintf("%08i", run_info.field_info.i_steps[i])).bp")
    # read the field
    file = adios_open_serial(fi, mode_readRandomAccess)
    for prop in properties
        field_name = "$(prop)_$(species)"
        if !(field_name in run_info.particle_info.properties)
            error("Property $prop for species $species not found in particle data.")
        end
        data[prop] = adios_load(file, "p$field_name", 0)
    end
    # remember to close the file
    close(file)

    return data
end