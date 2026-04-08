"""
    bp_to_vtk(data::EntityData, vtk_file::String, i_step::Integer)

Writes the fields of a given time step from an EntityData object to a VTK file.
"""
function bp_to_vtk(data::EntityData, vtk_file::String, i_step::Integer)

    dx = (data.settings["grid.extent"][2] - data.settings["grid.extent"][1]) / data.settings["grid.resolution"][1] + 1
    x = data.settings["grid.extent"][1]:dx:data.settings["grid.extent"][2]

    dy = (data.settings["grid.extent"][4] - data.settings["grid.extent"][3]) / data.settings["grid.resolution"][2] + 1
    y = data.settings["grid.extent"][3]:dy:data.settings["grid.extent"][4]

    dz = (data.settings["grid.extent"][6] - data.settings["grid.extent"][5]) / data.settings["grid.resolution"][3] + 1
    z = data.settings["grid.extent"][5]:dz:data.settings["grid.extent"][6]

    fields = read_field(data, data.fields.fields, i_step = i_step, verbose=true)

    vtk_grid(vtk_file, x, y, z) do vtk
        for field in data.fields.fields
            println("writing $field")
            vtk[field] = fields[field]
        end
    end
end
