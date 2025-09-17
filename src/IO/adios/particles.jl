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

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    for prop in properties
        field_name = "$(prop)_$(species)"
        if !(field_name in ["$(pr)_$(sp)" for pr in run_info.particles.properties, sp in run_info.particles.species])
            error("Property $prop for species $species not found in particle data.")
        end
        data[prop] = Vector{Float64}(undef, length(i))
    end 

    for j in i
        # construct the file name
        fi = joinpath(run_info.path, "particles", "particles.$(@sprintf("%08i", run_info.steps[j])).bp")
        # read the field
        file = adios_open_serial(fi, mode_readRandomAccess)
        for prop in properties
            field_name = "$(prop)_$(species)"
            data[prop][j] = reduction_function(adios_load(file, "p$field_name", 0))
        end    # remember to close the file
        close(file)

    end

    return data
end

"""
    read_particles( run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose=true )

Read `properties` for particle `species`. 
Give either file number `i`, time `t` or step number `i_step`.
If a `reduction_function` is given it computes the reduction per timestep.
In this case you can also give `i`, `t` or `i_step` as arrays to read time series.
"""
function read_particles(run_info::EntityData, species::Integer, 
                        properties::Vector{String}; 
                        i_step=nothing, t=nothing, i=nothing,
                        reduction_function=nothing,
                        verbose=true)

    if isnothing(i_step) && isnothing(t) && isnothing(i)
        error("Must provide either i_step, time t, or index i to read a field.")
    end

    if !isnothing(reduction_function)
        return read_particles_reduced(run_info, species, properties; 
                    i_step=i_step, t=t, i=i,
                    reduction_function=reduction_function,
                    verbose=verbose)
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

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    # construct the file name
    fi = joinpath(run_info.path, "particles", "particles.$(@sprintf("%08i", run_info.steps[i])).bp")
    # read the field
    file = adios_open_serial(fi, mode_readRandomAccess)
    for prop in properties
        field_name = "$(prop)_$(species)"
        if !(field_name in ["$(pr)_$(sp)" for pr in run_info.particles.properties, sp in run_info.particles.species])
            error("Property $prop for species $species not found in particle data.")
        end
        data[prop] = adios_load(file, "p$field_name", 0)
    end
    # remember to close the file
    close(file)

    return data
end