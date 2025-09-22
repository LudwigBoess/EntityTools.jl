"""
    read_spectrum(run_info::EntityData, species::Integer; 
                  i_step=nothing, t=nothing, i=nothing,
                  verbose=true)

Read spectrum for particle `species`. 
Give either file number `i`, time `t` or step number `i_step`.
"""
function read_spectrum(run_info::EntityData, species::Integer; 
                        i_step=nothing, t=nothing, i=nothing,
                        verbose=true)

    # get the relevant index/indeces for the file
    i = find_idxs(run_info, i_step, t, i, nothing, verbose)

    data = Dict()
    data["t"] = run_info.times[i]
    data["step"] = run_info.steps[i]

    # construct the file name
    fi = joinpath(run_info.path, "spectrum", "spectrum.$(@sprintf("%08i", run_info.steps[i])).bp")
    # read the field
    file = adios_open_serial(fi, mode_readRandomAccess)
    data["E"] = adios_load(file, "sEbn", 0)
    data["N"] = adios_load(file, "sN_$species", 0)
    # remember to close the file
    close(file)

    return data
end