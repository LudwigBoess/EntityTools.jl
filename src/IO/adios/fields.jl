function read_field(run_info::EntityData, field_name::String; 
                    i_step=nothing, t=nothing, i=nothing)

    if isnothing(i_step) && isnothing(t) && isnothing(i)
        error("Must provide either i_step, time t, or index i to read a field.")
    end

    if !isnothing(i_step)
        i = findfirst(isequal(i_step), run_info.field_info.i_steps)
        if isnothing(i)
            error("Step $i_step not found in field data.")
        end
    end
    if !isnothing(t)
        i = argmin(abs.(run_info.field_info.times .- t))
    end

    if i < 1 || i > length(run_info.field_info.i_steps)
        error("Index i=$i is out of bounds. Must be between 1 and $(length(run_info.field_info.i_steps)).")
    end

    # construct the file name
    fi = joinpath(RunInfo.path, "fields", "fields.$(@sprintf("%08i", run_info.field_info.i_steps[i])).bp")
    # read the field
    file = adios_open_serial(fi, mode_readRandomAccess)
    field = adios_load(file, "f$field_name", 0)
    # remember to close the file
    close(file)

    return field
end