"""
    read_fields_reduced(run_info::EntityData, field_names::Vector{String}; 
                    i, reduction_function, slice, verbose)

Reads fields over an output range and computes a reduction based on `reduction_function`.
"""
function read_fields_reduced(run_info::EntityData, field_names::Vector{String}; 
                    i, reduction_function, slice, verbose)

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    for field_name in field_names
        if !(field_name in run_info.fields.fields)
            error("Field $field_name not found in data.")
        end
        data[field_name] = Vector{Float64}(undef, length(i))
    end 

    # loop over files
    for (j, _step) in enumerate(i)
        # read the data of the current timestep
        step_data = read_field(run_info, field_names; 
                    i_step=run_info.steps[_step],
                    slice=slice,
                    verbose=verbose)

        # reduce the fields
        for field_name in field_names
            println("Reducing field $field_name at step $(j)...")
            data[field_name][j] = reduction_function(step_data[field_name])
        end
    end

    return data
end


"""
    read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    slice=nothing,
                    verbose::Bool=false)

Read fields with given `field_names` at a requested step `i_step`, time `t` or index `i`.
If a `reduction_function` it is applied to the fields in the requested step, time or index range.
"""
function read_field(run_info::EntityData, field_names::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    slice=nothing,
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
                    i, reduction_function, slice, verbose)
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

    # define selection
    start_arr = Vector{Int}(undef, length(Ndims))
    count_arr = Vector{Int}(undef, length(Ndims))

    for (i, Ndim) in enumerate(Ndims)
        X = adios_load(file, "X$(Ndim)", 0)
        if isnothing(slice)
            start_arr[i] = firstindex(X)-1 # ADIOS uses 0-based indexing
            count_arr[i] = length(X)
        else
            start_arr[i] = findfirst(X .>= slice[i][1])-1
            count_arr[i] = findlast(X .<= slice[i][2]) - start_arr[i]
        end
    end

    start_coords = Tuple(start_arr) 
    count_dims   = Tuple(count_arr)

    # read grid positions
    for (i, Ndim) in enumerate(Ndims)
        # select the variable to read
        grid_var = inquire_variable(file.io, "X$(Ndim)")
        # prepare selection
        set_selection(grid_var, start_coords, count_dims)
        # Dynamic allocation based on variable type
        T = type(grid_var)
        data["X$(Ndim)"] = Array{T}(undef, count_dims[i])
        get(file.engine, grid_var, data["X$(Ndim)"])
    end

    # read requested fields
    for field_name in field_names
        grid_var = inquire_variable(file.io, "f$field_name")
        # prepare selection
        set_selection(grid_var, start_coords, count_dims)
        # Dynamic allocation based on variable type
        T = type(grid_var)
        data["$field_name"] = Array{T}(undef, count_dims...)
        get(file.engine, grid_var, data["$field_name"])
    end

    # perform the reads
    perform_gets(file.engine)

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
            slice=nothing,
            verbose::Bool=false) = read_field(run_info, [field_name];
            i_step, t, i, reduction_function, slice, verbose)