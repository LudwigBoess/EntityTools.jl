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
                    slice=nothing,
                    verbose=true)

    

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    for prop in properties
        data[prop] = Vector{Float64}(undef, length(i))
    end 

    # loop over files
    for (j, _step) in enumerate(i)
        # read the data of the current timestep
        step_data = read_particles(run_info, species, properties; 
                    i_step=run_info.steps[_step],
                    slice=slice,
                    verbose=verbose)

        # reduce the particles
        for prop in properties
            data[prop][j] = reduction_function(step_data[prop])
        end    # remember to close the file
    end

    return data
end


"""
    reduce_read_positions(sel::Array{<:Integer})

Reduces the individual read positions by finding sub-ranges. 
Returns an array of read indices and an array with the number of entries per index.
"""
function reduce_read_positions(sel)

    Npart     = size(sel,1)
    index     = zeros(Int64, Npart) 
    n_to_read = zeros(Int64, Npart)

    index[1]     = sel[1] - 1
    n_to_read[1] = 1
    Nentries     = 1

    @inbounds for  i=2:Npart
        if (sel[i] == index[Nentries] + n_to_read[Nentries] + 1)
            n_to_read[Nentries] += 1
        else
            Nentries            += 1
            index[Nentries]     = sel[i] - 1
            n_to_read[Nentries] += 1
        end
    end

    resize!(index, Nentries)
    resize!(n_to_read, Nentries)

    return index, n_to_read
end


"""
    read_particles( run_info::EntityData, species::Integer, 
                    properties::Vector{String}; 
                    i_step=nothing, t=nothing, i=nothing,
                    reduction_function=nothing,
                    slice=nothing,
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
                        slice=nothing,
                        verbose=true)

    for prop in properties
        field_name = "$(prop)_$(species)"
        if !(field_name in ["$(pr)_$(sp)" for pr in run_info.particles.properties, sp in run_info.particles.species])
            error("Property $prop for species $species not found in particle data.")
        end
    end

    # get the relevant index/indeces for the file
    i = find_idxs(run_info, i_step, t, i, reduction_function, verbose)

    if !isnothing(reduction_function)
        return read_particles_reduced(run_info, species, properties; 
                    i=i,
                    reduction_function=reduction_function,
                    slice=slice,
                    verbose=verbose)
    end

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    # construct the file name
    fi = joinpath(run_info.path, "particles", "particles.$(@sprintf("%08i", run_info.steps[i])).bp")
    # open the file
    file = adios_open_serial(fi, mode_readRandomAccess)

    variables = adios_all_variable_names(file)
    X_fields = variables[startswith.(variables, "pX")]
    Ndims = unique([X[3] for X in X_fields])

    if isnothing(slice)
        slice = [ [-Inf, Inf] for _ in 1:length(Ndims) ]
    end

    # define selection
    selections = Vector{Vector{Int}}(undef, length(Ndims))
    for (i, Ndim) in enumerate(Ndims)
        if slice[i] == [-Inf, Inf]
            # no selection, take all particles
            part_var = inquire_variable(file.io, "pX$(Ndim)_$(species)")
            selections[i] = 1:Int64(ADIOS2.shape(part_var)[1])
        else
            # load the position array to find selected particles
            X = adios_load(file, "pX$(Ndim)_$(species)", 0)
            selections[i] = findall(slice[i][1] .<= X .<= slice[i][2])
        end
    end

    # remove unfiltered dimensions
    selections = selections[findall(slice != [-Inf, Inf])]

    # check which particles pass all selections
    sel = nothing
    if length(selections) == 1
        sel = selections[1]
    elseif length(Ndims) == 2
        sel = intersect(selections[1], selections[2])
    elseif length(Ndims) == 3
        sel = intersect(selections[1], selections[2], selections[3])
    end

    if verbose
        println("selected $(length(sel)) particles after applying spatial filters.")
    end
    # construct start positions and counts for bulk read
    start_arr, count_arr = reduce_read_positions(sel)
    # total number of particles to read
    n_to_read = sum(count_arr)

    start_coords = Tuple(start_arr) 
    count_dims   = Tuple(count_arr)

    for prop in properties
        field_name = "$(prop)_$(species)"
        if verbose
            println("preparing property $prop for species $species...")
        end
        part_var = inquire_variable(file.io, "p$field_name")
        # allocate array
        T = type(part_var)
        data[prop] = Array{T}(undef, n_to_read)

        N_read = 0
        for j in 1:length(count_arr)
            # prepare selection
            set_selection(part_var, (start_arr[j],), (count_arr[j],))
            # prepare the read
            get(file.engine, part_var, view(data[prop], N_read+1:N_read+count_arr[j]))
            N_read += count_arr[j]
        end
    end


    if verbose
        println("performing reads...")
    end
    # perform the reads
    perform_gets(file.engine)

    # remember to close the file
    close(file)

    if verbose
        println("done reading particle data.")
    end

    return data
end