# Disk-space estimators for Entity ADIOS2 output.
#
# These mirror exactly what the Entity C++ writers put on disk:
#   - checkpoints : src/framework/containers/fields_io.cpp  +  particles_io.cpp
#                   (declared/written in metadomain_chckpt.cpp)
#   - field    snapshots : src/output/writer.cpp  (Writer::defineFieldOutputs / WriteField)
#   - particle snapshots : src/framework/containers/particles_io.cpp (Particles::OutputWrite)
#
# Fundamental type sizes (src/global/global.h):
#   real_t   = double (8 B) or float (4 B) under -DSINGLE_PRECISION
#   prtldx_t = float  (4 B)  -- ALWAYS single precision, regardless of real_t
#   npart_t  = uint32_t (4 B)   ncells_t = uint32_t (4 B)
#   int      = 4 B           short (particle tag) = 2 B

# byte sizes of the fixed-width Entity types
const _BYTES_INT    = 4   # int            (cell index i, i_prev)
const _BYTES_SHORT  = 2   # short          (particle tag)
const _BYTES_PRTLDX = 4   # prtldx_t/float (sub-cell offset dx, dx_prev)
const _BYTES_NPART  = 4   # npart_t/uint32 (pld_i, tracking IDs)

"""
size of `real_t` in bytes for the requested `precision`. Accepts a Julia type
(`Float32`/`Float64`) or a symbol (`:float`/`:single` → 4 B, `:double` → 8 B).
"""
function _real_bytes(precision)
    precision in (Float32, :float, :single) && return 4
    precision in (Float64, :double)         && return 8
    throw(ArgumentError("precision must be Float32/:float/:single or Float64/:double, got $precision"))
end

"convert bytes to TB. `binary=true` (default) → TiB (1024^4); `binary=false` → decimal TB (1e12)."
_bytes_to_TB(bytes::Real; binary::Bool = true) = bytes / (binary ? 1024.0^4 : 1.0e12)


"""
    checkpoint_size(n_cells, n_particles; kwargs...) -> Float64

Required disk space, **in TB**, for a single Entity checkpoint (`step-XXXXXXXX.bp`)
holding `n_cells` field cells and `n_particles` particles.

Mirrors `Fields::CheckpointWrite` (fields) and `Particles::CheckpointWrite` (particles).

# Arguments
- `n_cells`: number of field cells stored. The checkpoint stores fields *including*
  ghost zones, so for an exact figure pass the ghosted count. Use
  [`cells_with_ghosts`](@ref) to inflate an active-cell count.
- `n_particles`: total number of particles stored (summed over all species).

# Keyword arguments
- `dim = 3`: spatial dimension (1, 2 or 3). Each dimension adds `i, dx, i_prev, dx_prev`.
- `engine = :SRPIC`: `:SRPIC` stores `em` (6 comp); `:GRPIC` additionally stores
  `em0` (6 comp) and `cur` (3 comp) → 15 components per cell.
- `precision = Float32`: `real_t` precision — `Float32`/`:float`/`:single` (4 B) or
  `Float64`/`:double` (8 B). Matches Entity's `-DSINGLE_PRECISION` build flag.
- `cartesian = true`: set `false` for a 2D non-Cartesian (q/spherical) run, which
  additionally stores a `phi` coordinate per particle.
- `npld_r = 0`: number of real-valued particle payloads per particle (`pld_r`).
- `npld_i = 0`: number of integer particle payloads per particle (`pld_i`, incl. tracking IDs).
- `binary = true`: report TiB (1024^4). Set `false` for decimal TB (1e12).

# Example
```julia
# 512^3 active cells (N_GHOSTS=2, single domain), 2e9 particles, double precision
nc = cells_with_ghosts((512, 512, 512))
checkpoint_size(nc, 2_000_000_000; dim = 3)
```
"""
function checkpoint_size(n_cells, n_particles;
                         dim::Integer    = 3,
                         engine::Symbol  = :SRPIC,
                         precision       = Float32,
                         cartesian::Bool = true,
                         npld_r::Integer = 0,
                         npld_i::Integer = 0,
                         binary::Bool    = true)

    (1 <= dim <= 3)             || throw(ArgumentError("dim must be 1, 2 or 3, got $dim"))
    engine in (:SRPIC, :GRPIC)  || throw(ArgumentError("engine must be :SRPIC or :GRPIC, got $engine"))

    R = _real_bytes(precision)

    # --- fields -------------------------------------------------------------
    # em (6 components) for SRPIC; + em0 (6) + cur (3) for GRPIC
    n_field_comp = engine === :GRPIC ? 15 : 6
    field_bytes  = float(n_cells) * n_field_comp * R

    # --- particles ----------------------------------------------------------
    # per dimension: i (int) + dx (prtldx) + i_prev (int) + dx_prev (prtldx)
    per_dim   = _BYTES_INT + _BYTES_PRTLDX + _BYTES_INT + _BYTES_PRTLDX  # 16 B
    # ux1, ux2, ux3 (3 reals) + weight (1 real)
    per_part  = dim * per_dim + 4 * R + _BYTES_SHORT                     # + tag (short)
    # phi coordinate only stored for 2D non-Cartesian
    (dim == 2 && !cartesian) && (per_part += R)
    # payloads
    per_part += npld_r * R + npld_i * _BYTES_NPART

    particle_bytes = float(n_particles) * per_part

    return _bytes_to_TB(field_bytes + particle_bytes; binary = binary)
end


# map a single T-tensor component character to Entity's integer code
# (out::InterpretComponents): t→0, x→1, y→2, z→3, digits 0..3 literally, i→-1, j→-2
function _comp_char_to_int(c::AbstractChar)
    c == 't' && return 0
    c == 'x' && return 1
    c == 'y' && return 2
    c == 'z' && return 3
    c == 'i' && return -1
    c == 'j' && return -2
    c in '0':'3' && return parse(Int, c)
    throw(ArgumentError("invalid T-tensor component '$c' (use t/x/y/z, 0..3, i or j)"))
end

# number of scalar components for a T-tensor name, replicating
# out::InterpretComponents with two indices (c1, c2)
function _T_component_count(c1::Int, c2::Int)
    if c1 < 0
        n = 0
        for i in 0:2
            if c2 < 0
                c2 == c1 ? (n += 1) : (n += length(i:2))  # diagonal vs sym. upper-tri
            else
                n += 1
            end
        end
        return n
    else
        return c2 < 0 ? 3 : 1
    end
end

"""
    field_components(field; engine = :SRPIC) -> Int

Number of scalar components Entity writes for a field-output `field`, replicating
`out::OutputField` (src/output/fields.cpp). `field` is a single name (e.g. `"E"`,
`"Tij"`, `"Rho_1_2"`) or a collection of names; for a collection the per-field counts
are summed.

Mapping:
- `E, B, D, H, J` → 3   (all three vector components)
- `V` → 3 for `:SRPIC`, 4 for `:GRPIC`
- `A` → 1   (only `A₃`; GRPIC only)
- `T…` → 1–6 from the index spec: `Ttt`/`T02`/`Txy` → 1, `T0i`/`Ti0` → 3, `Tij` → 6
- `Rho, Charge, N, Nppc, divE, divD`, and any custom name → 1

A trailing species selector (`_1`, `_1_2`, …) does not change the component count.
"""
function field_components(field::AbstractString; engine::Symbol = :SRPIC)
    engine in (:SRPIC, :GRPIC) ||
        throw(ArgumentError("engine must be :SRPIC or :GRPIC, got $engine"))

    # strip species suffix, then resolve the base field id (cf. OutputField ctor)
    base = String(split(field, '_')[1])
    lb   = lowercase(base)
    if lb == "dive" || lb == "divd"
        id = lb
    else
        # base id = everything before the first component/index character
        cut = findfirst(c -> c in ('0', '1', '2', '3', 'i', 'j', 'x', 'y', 'z', 't'),
                        base)
        id = lowercase(cut === nothing ? base : base[1:prevind(base, cut)])
    end

    if id in ("e", "b", "d", "h", "j")            # vector field / current
        return 3
    elseif id == "v"                              # bulk velocity
        return engine === :GRPIC ? 4 : 3
    elseif id == "a"                              # vector potential A3
        return 1
    elseif id == "t"                              # energy-momentum tensor
        length(field) >= 3 || throw(ArgumentError(
            "T-tensor field \"$field\" needs two index characters, e.g. \"Tij\", \"T0i\", \"Ttt\""))
        c1 = _comp_char_to_int(lowercase(field[2]))
        c2 = _comp_char_to_int(lowercase(field[3]))
        return _T_component_count(c1, c2)
    else                                          # scalars: Rho, divE/divD, N, Nppc, Charge, Custom
        return 1
    end
end

field_components(fields; engine::Symbol = :SRPIC) =
    sum(field_components(String(f); engine = engine) for f in fields)


"""
    field_snapshot_size(n_cells, fields; kwargs...) -> Float64

Required disk space, **in TB**, for one field snapshot (`fields/fields.XXXXXXXX.bp`)
over `n_cells` (full-resolution) active cells, storing the field outputs named in
`fields` — exactly as listed under `output.fields.quantities` in the Entity input,
e.g. `["E", "B", "J", "divE"]`.

Mirrors `Writer::defineFieldOutputs` / `WriteField`: each scalar component is written
as one `real_t` per output cell. The component count of each name is resolved with
[`field_components`](@ref) (`E/B/J` → 3, `divE` → 1, `Tij` → 6, …). Downsampling
(`output.fields.downsampling`) keeps every `downsample`-th cell per axis, reducing
the written cell count by `downsample^dim`.

# Arguments
- `n_cells`: number of full-resolution active cells (before downsampling). With
  `output.fields.ghosts = true` pass the ghosted count instead.
- `fields`: a field name, a collection of names (`["E","B","J","divE"]`), or an
  integer total number of scalar components.

# Keyword arguments
- `engine = :SRPIC`: affects only `V` (3 components for `:SRPIC`, 4 for `:GRPIC`).
- `downsample = 1`: downsampling factor. A scalar is applied per axis (effective
  reduction `downsample^dim`); a per-dimension tuple/vector reduces by `prod(downsample)`.
- `dim = 3`: spatial dimension, used only for a scalar `downsample`.
- `precision = Float32`: `real_t` precision — `Float32`/`:float`/`:single` (4 B) or
  `Float64`/`:double` (8 B).
- `binary = true`: report TiB (1024^4); `false` → decimal TB (1e12).

# Example
```julia
field_snapshot_size(1024^3, ["E", "B", "J", "divE"])           # 10 components
field_snapshot_size(1024^3, ["E", "B", "Tij"]; downsample = 2) # 3+3+6 = 12 components
```

Mesh coordinate arrays (`X1`, `X1e`, …) are negligible for production grids and
are not included.
"""
function field_snapshot_size(n_cells, fields;
                             engine::Symbol = :SRPIC,
                             downsample     = 1,
                             dim::Integer   = 3,
                             precision      = Float32,
                             binary::Bool   = true)
    (1 <= dim <= 3) || throw(ArgumentError("dim must be 1, 2 or 3, got $dim"))
    n_comp = fields isa Integer ? fields : field_components(fields; engine = engine)
    R      = _real_bytes(precision)
    dwn    = downsample isa Number ? float(downsample)^dim : prod(float.(downsample))
    n_out  = float(n_cells) / dwn
    return _bytes_to_TB(n_out * n_comp * R; binary = binary)
end


"""
    particle_snapshot_size(n_particles, stride; kwargs...) -> Float64

Required disk space, **in TB**, for one particle snapshot
(`particles/particles.XXXXXXXX.bp`) given `n_particles` particles written with output
`stride` (`output.particles.stride`).

Mirrors `Particles::OutputWrite`, which writes `nout = n_particles ÷ stride` particles,
each with weight `pW`, momenta `pU1..pU3`, and one position coordinate per dimension
(`pX1..`), all as `real_t`. Optional real payloads (`pPLDR*`) are `real_t`; tracking
IDs (`pIDX`, and `pRNK` under MPI) and integer payloads (`pPLDI*`) are `npart_t` (4 B).

# Arguments
- `n_particles`: number of particles before striding (per species, summed if you pass a total).
- `stride`: output stride; only every `stride`-th particle is written.

# Keyword arguments
- `dim = 3`: spatial dimension (1, 2 or 3).
- `precision = Float32`: `real_t` precision — `Float32`/`:float`/`:single` (4 B) or
  `Float64`/`:double` (8 B).
- `cartesian = true`: set `false` for a 2D non-Cartesian run (extra `pX3` phi coordinate).
- `npld_r = 0`: number of real payloads written per particle.
- `tracking = false`: whether species tracking is on (writes `pIDX`, plus `pRNK` under MPI).
- `mpi = true`: under MPI a tracked species also writes `pRNK`. Ignored if `tracking = false`.
- `npld_i = 0`: number of *additional* integer payloads written (beyond the tracking IDs).
- `binary = true`: report TiB (1024^4); `false` → decimal TB (1e12).
"""
function particle_snapshot_size(n_particles, stride;
                                dim::Integer    = 3,
                                precision       = Float32,
                                cartesian::Bool = true,
                                npld_r::Integer = 0,
                                tracking::Bool  = false,
                                mpi::Bool       = true,
                                npld_i::Integer = 0,
                                binary::Bool    = true)

    (1 <= dim <= 3) || throw(ArgumentError("dim must be 1, 2 or 3, got $dim"))
    stride >= 1     || throw(ArgumentError("stride must be >= 1, got $stride"))

    R    = _real_bytes(precision)
    nout = floor(float(n_particles) / stride)

    # number of position coordinates written: one per dim, +1 (phi) for 2D non-Cartesian
    n_pos = dim + ((dim == 2 && !cartesian) ? 1 : 0)

    # real_t quantities: weight (1) + momenta (3) + positions + real payloads
    per_part = (1 + 3 + n_pos + npld_r) * R
    # tracking IDs (npart_t)
    if tracking
        per_part += _BYTES_NPART                 # pIDX
        mpi && (per_part += _BYTES_NPART)         # pRNK
    end
    # additional integer payloads (npart_t)
    per_part += npld_i * _BYTES_NPART

    return _bytes_to_TB(nout * per_part; binary = binary)
end


"""
    cells_with_ghosts(n_active; n_ghosts = 2, n_domains = 1) -> Int

Inflate a per-dimension active-cell count to the ghosted cell count that a checkpoint
actually stores (see `Metadomain::InitCheckpointWriter`):

    glob_shape_with_ghosts[d] = n_active[d] + 2 * n_ghosts * n_domains_per_dim[d]

and returns the product over dimensions (total stored cells).

# Arguments
- `n_active`: integer, or a tuple/vector of per-dimension active cell counts.
- `n_ghosts = 2`: `N_GHOSTS` (2 for `SHAPE_ORDER = 2`).
- `n_domains = 1`: integer or per-dimension tuple/vector of domain (MPI) counts per axis.

# Example
```julia
cells_with_ghosts((512, 512, 512))            # single domain, N_GHOSTS=2
cells_with_ghosts((512, 512, 512); n_domains = (4, 4, 2))
```
"""
function cells_with_ghosts(n_active; n_ghosts::Integer = 2, n_domains = 1)
    nd = n_active isa Number ? (n_active,) : Tuple(n_active)
    dm = n_domains isa Number ? ntuple(_ -> n_domains, length(nd)) : Tuple(n_domains)
    length(dm) == length(nd) ||
        throw(ArgumentError("n_domains must match the number of dimensions in n_active"))
    total = 1
    for d in eachindex(nd)
        total *= nd[d] + 2 * n_ghosts * dm[d]
    end
    return total
end
