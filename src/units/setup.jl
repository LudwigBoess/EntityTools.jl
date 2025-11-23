"""
    EntityUnits(data::EntityData)

Constructs the units from the simulation settings stored in `data`.
Returns an [`EntityUnits`](@ref) struct containing the normalized units for the simulation.
"""
function EntityUnits(data::EntityData)
    
    return EntityUnits( data.settings["scales.skindepth0"], 
                        data.settings["scales.larmor0"],
                        data.settings["scales.sigma0"],
                        data.settings["scales.B0"],
                        data.settings["scales.omegaB0"],
                        data.settings["scales.q0"],
                        data.settings["scales.V0"],
                        data.settings["scales.n0"] )
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