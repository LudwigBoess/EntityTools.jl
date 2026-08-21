# =============================================================================
#  Reading Entity's [scales] out of a run
# =============================================================================
#
#  Everything the unit chain needs is stored in the ADIOS2 attributes of an
#  output file, with two exceptions:
#
#    - the species masses (`particles.species` is not written as an attribute),
#      which are recovered from the `<name>.info` report Entity writes next to
#      the output folders, or read straight from the parameter file;
#    - the injected density fraction `n_frac`, which lives in the problem
#      generator and has to be supplied (see [`EntityUnits`](@ref)).
# =============================================================================

"""
    read_species(info_path::AbstractString)

Species table parsed from the `<name>.info` report Entity writes into the run
directory. Returns a `Vector` of `(; label, mass, charge)` in the order the
species are declared, or an empty vector if the file is missing or holds no
species block.

The masses are the only place they can be recovered from a finished run:
`particles.species` is not written to the output attributes.
"""
function read_species(info_path::AbstractString)

    species = NamedTuple{(:label, :mass, :charge), Tuple{String, Float64, Float64}}[]
    isfile(info_path) || return species

    lines = readlines(info_path)
    # a species block opens with a bare `- Species #N:`; the per-domain memory
    # footprint uses `- Species #N (label): <size>` and must not match
    header = r"^\s*-\s*Species\s+#\d+:\s*$"
    any_species = r"^\s*-\s*Species\s+#"

    i = 1
    while i ≤ length(lines)
        if occursin(header, lines[i])
            label, mass, charge = "", NaN, NaN
            j = i + 1
            while j ≤ length(lines) && !occursin(any_species, lines[j])
                m = match(r"^\s*-\s*Label:\s*(\S+)", lines[j])
                m !== nothing && (label = String(m[1]))
                m = match(r"^\s*-\s*Mass:\s*([-+.\deE]+)", lines[j])
                m !== nothing && (mass = parse(Float64, m[1]))
                m = match(r"^\s*-\s*Charge:\s*([-+.\deE]+)", lines[j])
                m !== nothing && (charge = parse(Float64, m[1]))
                j += 1
            end
            isnan(mass) || push!(species, (; label, mass, charge))
            i = j
        else
            i += 1
        end
    end
    return species
end

"""
    read_species(data::EntityData)

Species table of a run, read from `<path>/<simulation.name>.info`, falling back
to the only `.info` file in the run directory if that name does not exist.
"""
function read_species(data::EntityData)
    named = joinpath(data.path, string(get(data.settings, "simulation.name", "")) * ".info")
    isfile(named) && return read_species(named)
    infos = filter(f -> endswith(f, ".info"), readdir(data.path))
    length(infos) == 1 || return read_species("")
    return read_species(joinpath(data.path, only(infos)))
end

"""
    mass_ratio_from_species(species) -> Union{Float64, Nothing}

Ion/electron mass ratio μ = m_i/m_e of a species table: the heaviest over the
lightest charged species. Returns `1.0` for a single species (pair plasma
convention) and `nothing` if the table holds no charged species.
"""
function mass_ratio_from_species(species)
    charged = filter(s -> s.charge != 0 && s.mass > 0, species)
    isempty(charged) && return nothing
    return maximum(s.mass for s in charged) / minimum(s.mass for s in charged)
end

"""
    cell_size(settings, dim, V0)

Fiducial cell size [code length] of a run. Taken from the grid extent and
resolution for a Cartesian grid, and from the cell volume otherwise.

Entity's own `scales.dx0` attribute is *not* the cell size: it is `dxMin`, the
CFL length, which in Minkowski is dx/sqrt(dim).
"""
function cell_size(settings::AbstractDict, dim::Integer, V0::Real)

    cart = string(get(settings, "grid.metric.coord", "cart")) == "cart"
    if cart && haskey(settings, "grid.extent") && haskey(settings, "grid.resolution")
        # `grid.extent` is stored flat as [x1min, x1max, x2min, x2max, ...]
        extent = Float64.(settings["grid.extent"])
        res    = Float64.(settings["grid.resolution"])
        dxs    = [(extent[2d] - extent[2d-1]) / res[d] for d in 1:dim]
        # non-cubic cells make "the" cell size ambiguous, so say so and fall back
        if all(isapprox(dxs[1]; rtol = 1e-4), dxs)
            return dxs[1]
        end
        @warn "non-cubic cells $(dxs): using the geometric mean V0^(1/dim) as dx"
    elseif !cart
        @warn "non-Cartesian grid: dx is the fiducial cell volume^(1/dim), " *
              "and all lengths are code coordinates"
    end
    return Float64(V0)^(1 / dim)
end

"""
    EntityUnits(data::EntityData; kwargs...)

Construct the unit conversion factors from the settings stored in the output of
a run.

`skindepth0`, `larmor0`, `V0`, `ppc0` and the grid geometry come from the ADIOS2
attributes. The mass ratio is taken from the `<name>.info` report and the
temperature from `setup.temperature`/`setup.temperature_ratio` when the problem
generator defines them; any of them can be overridden by keyword.

`n_frac` cannot be recovered from the output — it is set by the problem
generator's injection call. It defaults to `1` (the fiducial scales Entity
reports); the standard two-species SRPIC injection needs `n_frac = 0.5`.

# Keyword arguments
- `mass_ratio = nothing`        : μ = m_i/m_e. `nothing` reads it from `<name>.info`.
- `n_frac = 1`                  : n_s/n0 for one species.
- `temperature = nothing`       : T_e [m_e c²]. `nothing` reads `setup.temperature`.
- `temperature_ratio = nothing` : T_i/T_e. `nothing` reads `setup.temperature_ratio`.
- `mass_density = :ion`         : mass density behind v_A and c_s.
- `gamma_ad = 5/3`              : adiabatic index behind c_s.
- `check = true`                : cross-check the derived `sigma0`, `n0` and `q0`
  against the values Entity stored, and warn on a mismatch.
"""
EntityUnits(data::EntityData; kwargs...) =
    _units_from_settings(data.settings, read_species(data); source = data.path, kwargs...)

"""
    check_scales(u::EntityUnits, settings::AbstractDict; rtol = 1e-4)

Compare the fiducial quantities re-derived by [`EntityUnits`](@ref) against the
ones Entity stored in the output attributes, and warn on any mismatch. Guards
against a change in Entity's own normalisation going unnoticed.
"""
function check_scales(u::EntityUnits, settings::AbstractDict; rtol::Real = 1e-4)
    for (key, derived) in (("scales.sigma0", u.sigma0), ("scales.B0", u.B0),
                           ("scales.omegaB0", u.omegaB0), ("scales.n0", u.n0),
                           ("scales.q0", u.q0))
        haskey(settings, key) || continue
        stored = Float64(settings[key])
        isapprox(stored, derived; rtol) ||
            @warn "$key = $stored in the output but $derived from the unit chain"
    end
    return nothing
end

"""
    EntityUnits(path::String; kwargs...)

Construct the unit conversion factors from a run directory or from a parameter
file.

A directory is opened as an [`EntityData`](@ref); a `.toml` file is parsed
directly, which also recovers the species masses (and hence the mass ratio) from
`[[particles.species]]`. Keywords are those of
[`EntityUnits(::EntityData)`](@ref).
"""
function EntityUnits(path::String; kwargs...)

    isdir(path) && return EntityUnits(EntityData(path); kwargs...)
    isfile(path) || error("$path is neither a directory nor a file!")

    data = TOML.parsefile(path)
    settings = Dict{String, Any}()

    # flatten the parts of the parameter file the unit chain needs into the same
    # key layout the output attributes use
    settings["scales.skindepth0"] = Float64(data["scales"]["skindepth0"])
    settings["scales.larmor0"]    = Float64(data["scales"]["larmor0"])
    settings["particles.ppc0"]    = Float64(data["particles"]["ppc0"])

    grid = data["grid"]
    resolution = Int.(grid["resolution"])
    dim        = length(resolution)
    # in the parameter file the extent is nested pairwise; store it flat
    extent = Float64[e for pair in grid["extent"] for e in pair]
    settings["grid.dim"]          = dim
    settings["grid.extent"]       = extent
    settings["grid.resolution"]   = resolution
    settings["grid.metric.coord"] = get(get(grid, "metric", Dict()), "metric", "minkowski") ==
                                    "minkowski" ? "cart" : "curvilinear"
    settings["scales.V0"] = prod((extent[2d] - extent[2d-1]) / resolution[d] for d in 1:dim)

    setup = get(data, "setup", Dict())
    haskey(setup, "temperature") &&
        (settings["setup.temperature"] = Float64(setup["temperature"]))
    haskey(setup, "temperature_ratio") &&
        (settings["setup.temperature_ratio"] = Float64(setup["temperature_ratio"]))

    species = [(; label = String(get(sp, "label", "")),
                  mass = Float64(get(sp, "mass", NaN)),
                  charge = Float64(get(sp, "charge", 0.0)))
               for sp in get(data["particles"], "species", [])]

    return _units_from_settings(settings, species; source = path, kwargs...)
end

# Shared back end of both constructors: the settings are in the flat key layout
# of the output attributes, the species table is supplied separately because it
# comes from a different place in each case.
function _units_from_settings(settings::AbstractDict, species;
                              mass_ratio           = nothing,
                              n_frac::Real         = 1,
                              temperature          = nothing,
                              temperature_ratio    = nothing,
                              mass_density::Symbol = :ion,
                              gamma_ad::Real       = 5 / 3,
                              check::Bool          = true,
                              source::AbstractString = "")

    d0  = Float64(settings["scales.skindepth0"])
    r0  = Float64(settings["scales.larmor0"])
    V0  = Float64(settings["scales.V0"])
    p0  = Float64(settings["particles.ppc0"])
    dim = Int(settings["grid.dim"])
    dx  = cell_size(settings, dim, V0)

    # μ: the species masses are not among the output attributes
    μ = mass_ratio === nothing ? mass_ratio_from_species(species) : Float64(mass_ratio)
    if μ === nothing
        Int(get(settings, "particles.nspec", length(species))) > 1 &&
            @warn "could not read the species masses" * (isempty(source) ? "" : " of $source") *
                  "; assuming mass_ratio = 1. Pass `mass_ratio` explicitly."
        μ = 1.0
    end

    # T_e / T_i: only present when the problem generator declares them
    Te = temperature === nothing ? Float64(get(settings, "setup.temperature", NaN)) :
                                   Float64(temperature)
    τ  = temperature_ratio === nothing ?
         Float64(get(settings, "setup.temperature_ratio", 1)) : Float64(temperature_ratio)

    u = EntityUnits(d0, r0; dim, dx, V0, ppc0 = p0, mass_ratio = μ, n_frac,
                    temperature = Te, temperature_ratio = τ, mass_density, gamma_ad)

    check && check_scales(u, settings)
    return u
end
