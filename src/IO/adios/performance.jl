"""
    adios2_throughput(path::AbstractString)

Read an ADIOS2 BP5 `profiling.json` and return a NamedTuple with:

- `aggregate_TBps`     : Σwbytes(all ranks) / max(durable_us) [decimal TB/s]
- `per_rank_mean_GBps` : mean over writer ranks of wbytes / durable_us [decimal GB/s]
- `per_rank_std_GBps`  : stddev of the same
- `per_rank_raw_GBps`  : mean of wbytes / Σwrite.mus over writer ranks [decimal GB/s] —
                         the raw POSIX write() rate, kept for reference; this measures
                         the rate of memcpy into the page/client cache and overstates
                         durable bandwidth, often above the filesystem spec
- `n_writers`          : number of ranks that issued writes (BP5 aggregators)
- `n_ranks`            : total ranks in the profile

`durable_us` per rank is `ES + ES_close + DC_WaitOnAsync1 + DC_WaitOnAsync2`,
i.e. the EndStep window plus the Close-time drain of BP5's async writer thread.
With `AsyncWrite=ON`, `EndStep` returns before bytes are durable; the drain shows
up as `DC_WaitOnAsync*` at Close. Using this window gives bandwidth bounded by
what actually reached storage.

Aggregate throughput uses `max(durable_us)` over all ranks (writers and
non-writers) — non-writer ranks still block at Close waiting for the aggregator,
and the wall clock is bounded by the slowest rank.
"""
function adios2_throughput(path::AbstractString)
    raw = read(path, String)
    # ADIOS2 sometimes emits a trailing comma before the closing array bracket
    raw = replace(raw, r",(\s*\])" => s"\1")
    data = JSON3.read(raw)

    total_bytes          = 0
    rank_eff_GBps        = Float64[]
    rank_raw_GBps        = Float64[]
    durable_times_us     = Float64[]

    for rec in data
        rank_bytes    = 0
        rank_write_us = 0
        for (k, v) in pairs(rec)
            startswith(String(k), "transport_") || continue
            v isa JSON3.Object || continue
            rank_bytes    += get(v, :wbytes, 0)
            wr             = get(v, :write, nothing)
            wr === nothing || (rank_write_us += get(wr, :mus, 0))
        end
        total_bytes += rank_bytes

        es      = Float64(get(rec, :ES_mus, 0))
        esclose = Float64(get(rec, :ES_close_mus, 0))
        wait1   = Float64(get(rec, :DC_WaitOnAsync1_mus, 0))
        wait2   = Float64(get(rec, :DC_WaitOnAsync2_mus, 0))
        durable_us = es + esclose + wait1 + wait2
        durable_us > 0 && push!(durable_times_us, durable_us)

        if rank_bytes > 0 && durable_us > 0
            # bytes / μs == MB/s (decimal); /1000 -> GB/s
            push!(rank_eff_GBps, rank_bytes / durable_us / 1000)
        end
        if rank_bytes > 0 && rank_write_us > 0
            push!(rank_raw_GBps, rank_bytes / rank_write_us / 1000)
        end
    end

    isempty(durable_times_us) && error("no ES/Close/Async timing entries found in $path")
    wall_us = maximum(durable_times_us)
    # bytes / μs == MB/s; /1e6 -> TB/s (decimal)
    aggregate_TBps = total_bytes / wall_us / 1e6

    return (
        aggregate_TBps     = aggregate_TBps,
        per_rank_mean_GBps = isempty(rank_eff_GBps) ? NaN : mean(rank_eff_GBps),
        per_rank_std_GBps  = length(rank_eff_GBps) < 2 ? NaN : std(rank_eff_GBps),
        per_rank_raw_GBps  = isempty(rank_raw_GBps) ? NaN : mean(rank_raw_GBps),
        n_writers          = length(rank_eff_GBps),
        n_ranks            = length(data),
    )
end
