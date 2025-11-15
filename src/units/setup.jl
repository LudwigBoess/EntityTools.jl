"""
    EntityUnits(data::EntityData)

Reads the normalized units from the first available file in the `fields`, `particles`, or `spectra` subfolders of the simulation path.
Returns an [`EntityUnits`](@ref) struct containing the normalized units for the simulation.
"""
function EntityUnits(data::EntityData)
    
    if !isnothing(data.fields)
        file = select_first_file(data.path, "fields")
    elseif !isnothing(data.particles)
        file = select_first_file(data.path, "particles")
    elseif !isnothing(data.spectra)
        file = select_first_file(data.path, "spectra")
    else
        error("No file found in fields, particles, or spectra subfolders.")
    end

    return EntityUnits( adios_attribute_data(file, "scales.skindepth0"), 
                        adios_attribute_data(file, "scales.larmor0"),
                        adios_attribute_data(file, "scales.sigma0"),
                        adios_attribute_data(file, "scales.B0"),
                        adios_attribute_data(file, "scales.omegaB0"),
                        adios_attribute_data(file, "scales.q0"),
                        adios_attribute_data(file, "scales.V0"),
                        adios_attribute_data(file, "scales.n0") )
end

"""
    EntityUnits(toml_path::String)

Reads the normalized units from the simulations' parameter file at the specified path.
"""
function EntityUnits(toml_path::String)
    if !isfile(toml_path)
        error("File $toml_path not found!")
    end
    data = TOML.parsefile(toml_path)
    extent = data["grid"]["extent"][1]

    return EntityUnits( data["scales"]["skindepth0"], 
                        data["scales"]["larmor0"],
                        data["particles"]["ppc0"],
                        data["grid"]["extent"],
                        data["grid"]["resolution"],
                        data["grid"]["metric"]["metric"] )
end