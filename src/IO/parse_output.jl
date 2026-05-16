# ---------------------------------------------------------------------------
# Time-unit handling
# ---------------------------------------------------------------------------

const TIME_UNITS = Dict(
    "ns"  => 1e-9,
    "µs"  => 1e-6,
    "us"  => 1e-6,
    "ms"  => 1e-3,
    "s"   => 1.0,
    "min" => 60.0,
    "h"   => 3600.0,
)

to_seconds(value::Real, unit::AbstractString) = value * TIME_UNITS[unit]

const TIME_PAIR = r"([\d.]+)\s*(ns|µs|us|ms|min|h|s)"

"""
    parse_compound_time(str)

Parse compound durations like `"10s 143ms"`, `"3min 39s"`, `"286ms 625µs"`
into seconds.
"""
function parse_compound_time(str::AbstractString)
    total = 0.0
    for m in eachmatch(TIME_PAIR, str)
        total += to_seconds(parse(Float64, m.captures[1]), m.captures[2])
    end
    return total
end

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

"""
    StepData

Per-step record from an entity `.out` file. All times are in seconds.
Substep timings are accessible by name: `step["CurrentDeposit"]`.
"""
struct StepData
    step::Int
    sim_time::Float64
    dt::Float64
    species::Dict{String,Float64}          # label => particle count (global total)
    n_active::Float64                      # Σ species counts
    substeps::Dict{String,Float64}         # substep name => seconds
    pre_total_substeps::Vector{String}     # order of substeps included in Total
    post_total_substeps::Vector{String}    # order of substeps reported after Total
    total::Float64                         # "Total" line in seconds
    timestep_duration::Float64             # seconds
    remaining_time::Float64                # seconds
    elapsed_time::Float64                  # seconds
end

Base.getindex(s::StepData, name::AbstractString) = s.substeps[name]
Base.haskey(s::StepData, name::AbstractString)   = haskey(s.substeps, name)

"""
    SimOutput

Container for all per-step records parsed from a `.out` file. Iterable and
indexable by step index. The set of substep and species names observed across
the whole run is kept in `substep_names` / `species_labels`.
"""
struct SimOutput
    steps::Vector{StepData}
    substep_names::Vector{String}          # union of substep names, insertion order
    species_labels::Vector{String}         # union of species labels, insertion order
    source::String                         # path the data was parsed from
end

Base.length(s::SimOutput)               = length(s.steps)
Base.getindex(s::SimOutput, i::Integer) = s.steps[i]
Base.iterate(s::SimOutput, st...)       = iterate(s.steps, st...)
Base.firstindex(s::SimOutput) = 1
Base.lastindex(s::SimOutput)  = length(s.steps)

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

const STEP_LINE   = r"^Step:\s*(\d+)"
const NUMBER      = raw"[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"
const TIME_LINE   = Regex("^Time:\\s*($NUMBER).*?Δt\\s*=\\s*($NUMBER)")
const SUBSTEP_LN  = r"^\s+([A-Za-z][A-Za-z0-9_]*)\.{2,}\s*([\d.]+)\s*(ns|µs|us|ms|min|h|s)\b"
const TOTAL_LINE  = r"^Total\s+([\d.]+)\s*(ns|µs|us|ms|min|h|s)\b"
const SPECIES_LN  = r"^\s+species\s+\d+\s+\(([^)]+)\)\.{2,}\s*([\d.eE+-]+)"
const TIMESTEPDUR = r"^Timestep duration:\s*(.+)$"
const REMAINING   = r"^Remaining time:\s*(.+)$"
const ELAPSED     = r"^Elapsed time:\s*(.+)$"

"""
    parse_output(path) -> SimOutput

Parse an entity `.out` log file. The number and names of substeps are
discovered dynamically. All times are stored in seconds.
"""
function parse_output(path::AbstractString)
    isfile(path) || error("file not found: $path")

    steps         = StepData[]
    substep_order = String[]
    species_order = String[]

    in_step       = false
    cur_step      = 0
    cur_simtime   = 0.0
    cur_dt        = 0.0
    cur_species   = Dict{String,Float64}()
    cur_substeps  = Dict{String,Float64}()
    cur_pre       = String[]
    cur_post      = String[]
    cur_total     = NaN
    cur_tsdur     = NaN
    cur_remaining = NaN
    cur_elapsed   = NaN
    saw_total     = false

    flush_step! = function ()
        n_active = isempty(cur_species) ? 0.0 : sum(values(cur_species))
        push!(steps, StepData(
            cur_step, cur_simtime, cur_dt,
            copy(cur_species), n_active,
            copy(cur_substeps), copy(cur_pre), copy(cur_post),
            cur_total, cur_tsdur, cur_remaining, cur_elapsed,
        ))
    end

    reset_step! = function ()
        empty!(cur_species); empty!(cur_substeps)
        empty!(cur_pre);     empty!(cur_post)
        cur_total = NaN; cur_tsdur = NaN
        cur_remaining = NaN; cur_elapsed = NaN
        saw_total = false
    end

    open(path, "r") do io
        for line in eachline(io)
            m = match(STEP_LINE, line)
            if m !== nothing
                in_step && flush_step!()
                reset_step!()
                in_step  = true
                cur_step = parse(Int, m.captures[1])
                continue
            end
            in_step || continue

            m = match(TIME_LINE, line)
            if m !== nothing
                cur_simtime = parse(Float64, m.captures[1])
                cur_dt      = parse(Float64, m.captures[2])
                continue
            end

            m = match(TOTAL_LINE, line)
            if m !== nothing
                cur_total = to_seconds(parse(Float64, m.captures[1]), m.captures[2])
                saw_total = true
                continue
            end

            m = match(SUBSTEP_LN, line)
            if m !== nothing
                name = m.captures[1]
                t    = to_seconds(parse(Float64, m.captures[2]), m.captures[3])
                cur_substeps[name] = t
                if saw_total
                    push!(cur_post, name)
                else
                    push!(cur_pre, name)
                end
                name in substep_order || push!(substep_order, name)
                continue
            end

            m = match(SPECIES_LN, line)
            if m !== nothing
                label = m.captures[1]
                cur_species[label] = parse(Float64, m.captures[2])
                label in species_order || push!(species_order, label)
                continue
            end

            m = match(TIMESTEPDUR, line)
            if m !== nothing
                cur_tsdur = parse_compound_time(m.captures[1])
                continue
            end

            m = match(REMAINING, line)
            if m !== nothing
                cur_remaining = parse_compound_time(m.captures[1])
                continue
            end

            m = match(ELAPSED, line)
            if m !== nothing
                cur_elapsed = parse_compound_time(m.captures[1])
                continue
            end
        end
        in_step && flush_step!()
    end

    return SimOutput(steps, substep_order, species_order, abspath(path))
end

# ---------------------------------------------------------------------------
# Series accessors
# ---------------------------------------------------------------------------

"""
    substep_series(sim, name)

Per-step time series (seconds) for substep `name`. Steps where the substep
is absent return `NaN`.
"""
substep_series(sim::SimOutput, name::AbstractString) =
    [get(s.substeps, name, NaN) for s in sim.steps]

total_series(sim::SimOutput)        = [s.total              for s in sim.steps]
active_particles(sim::SimOutput)    = [s.n_active           for s in sim.steps]
timestep_durations(sim::SimOutput)  = [s.timestep_duration  for s in sim.steps]
elapsed_times(sim::SimOutput)       = [s.elapsed_time       for s in sim.steps]
remaining_times(sim::SimOutput)     = [s.remaining_time     for s in sim.steps]
sim_times(sim::SimOutput)           = [s.sim_time           for s in sim.steps]
step_numbers(sim::SimOutput)        = [s.step               for s in sim.steps]

species_series(sim::SimOutput, label::AbstractString) =
    [get(s.species, label, NaN) for s in sim.steps]

# ---------------------------------------------------------------------------
# Aggregates
# ---------------------------------------------------------------------------

_clean(xs) = filter(!isnan, xs)
_mean(xs)  = isempty(xs) ? NaN : sum(xs) / length(xs)

"""
    sum(sim::SimOutput)              -> seconds
    sum(sim::SimOutput, name)        -> seconds
    maximum(sim::SimOutput, name)    -> seconds
    minimum(sim::SimOutput, name)    -> seconds
    mean(sim::SimOutput)             -> seconds
    mean(sim::SimOutput, name)       -> seconds

Aggregate over the simulation. With no `name`, the aggregate runs over the
per-step `Total` line. With a `name`, it runs over the named substep
(e.g. `"CurrentDeposit"`); steps where the substep is absent are skipped.
"""
Base.sum(sim::SimOutput) = sum(_clean(total_series(sim)); init = 0.0)

Base.sum(sim::SimOutput, name::AbstractString) =
    sum(_clean(substep_series(sim, name)); init = 0.0)

Base.maximum(sim::SimOutput, name::AbstractString) =
    maximum(_clean(substep_series(sim, name)); init = -Inf)

Base.minimum(sim::SimOutput, name::AbstractString) =
    minimum(_clean(substep_series(sim, name)); init =  Inf)

Statistics.mean(sim::SimOutput) = _mean(_clean(total_series(sim)))

Statistics.mean(sim::SimOutput, name::AbstractString) =
    _mean(_clean(substep_series(sim, name)))

"""
    summary_table(sim) -> Vector{NamedTuple}

One row per substep with `(sum, mean, min, max)` in seconds, in the order the
substeps first appeared in the file.
"""
function summary_table(sim::SimOutput)
    [(name = n,
      sum  = sum(sim, n),
      mean = Statistics.mean(sim, n),
      min  = minimum(sim, n),
      max  = maximum(sim, n)) for n in sim.substep_names]
end
