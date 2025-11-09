"""
    read_fields_reduced(run_info::EntityData, field_names::Vector{String}; 
                    i, reduction_function, verbose)

Reads fields over an output range and computes a reduction based on `reduction_function`.
"""
function read_fields_reduced(run_info::EntityData, field_names::Vector{String}; 
                    i, reduction_function, verbose)

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    for field_name in field_names
        if !(field_name in run_info.fields.fields)
            error("Field $field_name not found in data.")
        end
        data[field_name] = Vector{Float64}(undef, length(i))
    end 

    for j in i
        # construct the file name
        fi = joinpath(run_info.path, "fields", "fields.$(@sprintf("%08i", run_info.steps[j])).bp")
        # read the field
        file = adios_open_serial(fi, mode_readRandomAccess)
        for field_name in field_names
            data[field_name][j] = reduction_function(adios_load(file, "f$field_name", 0))
        end    
        # remember to close the file
        close(file)
    end

    return data
end


"""
    read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false)

Read fields with given `field_names` at a requested step `i_step`, time `t` or index `i`.
If a `reduction_function` it is applied to the fields in the requested step, time or index range.
"""
function read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    verbose::Bool=false)

    # check if all requested fields are present
    for field_name in field_names
        if !(field_name in run_info.fields.fields)
            error("Field $field_name not found in data.")
        end
    end

    # find relevant indices for requested times/steps
    i = find_idxs(run_info, i_step, t, i, reduction_function, verbose)

    # compute reduced fields over requested timesteps
    if !isnothing(reduction_function)
        return read_fields_reduced(run_info, field_names; 
                    i=i,
                    reduction_function=reduction_function,
                    verbose=verbose)
    end

    # allocate storage dict
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

Reads a single field for an individual file. 
If `reduction_function` is provied the reduced field is computed for the requested output or a range of outputs.
"""
read_field(run_info::EntityData, field_name::String; 
            i_step=nothing, t=nothing, i=nothing,
            reduction_function=nothing,
            verbose::Bool=false) = read_field(run_info, [field_name];
            i_step, t, i, reduction_function, verbose)