function read_particles_reduced(run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i, verbose)

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
    read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false)

Read fields with given `field_names`.
"""
function read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false)

    i = find_idxs(run_info, i_step, t, i, reduction_function, verbose)

    if !isnothing(reduction_function)
        return read_particles_reduced(run_info, species, properties; 
                    i=i,
                    reduction_function=reduction_function,
                    verbose=verbose)
    end

    for j = 1:Nfields
        if !(field_names[j] in run_info.fields.fields)
            error("Field $field_name not found in data.")
        end
    end

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]
    
    # construct the file name
    fi = joinpath(run_info.path, "fields", "fields.$(@sprintf("%08i", run_info.steps[i])).bp")
    # read the field
    file = adios_open_serial(fi, mode_readRandomAccess)
    variables = adios_all_variable_names(file)
    
    # read grid positions
    X_fields = variables[startswith.(variables, "X")]
    Ndims = unique([X[2] for X in X_fields])
    for Ndim in Ndims 
        data["X$(Ndim)e"] = adios_load(file, "X$(Ndim)e", 0)
    end

    # read requested fields
    for field_name in field_names
        data[field_name] = adios_load(file, "f$field_name", 0)
    end
    # remember to close the file
    close(file)

    return data
end


"""
    read_field(run_info::EntityData, field_name::String; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false)

Same
"""
read_field(run_info::EntityData, field_name::String; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false) = read_field(run_info, [field_name];
                    i_step, t, i, reduction_function, verbose)