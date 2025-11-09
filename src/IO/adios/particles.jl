"""
    read_particles_reduced(run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i=nothing, reduction_function=nothing,
                    verbose=true)

Read `properties` for particle `species` and apply `reduction_function` to the data.
"""
function read_particles_reduced(run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i=nothing, reduction_function=nothing,
                    verbose=true)

    

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

    # get the relevant index/indeces for the file
    i = find_idxs(run_info, i_step, t, i, reduction_function, verbose)

    if !isnothing(reduction_function)
        return read_particles_reduced(run_info, species, properties; 
                    i=i,
                    reduction_function=reduction_function,
                    verbose=verbose)
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