using Printf
using Unitful
using DelimitedFiles

"""
    reinterpret_time(time_string)

Takes the time string and computes the duration in seconds
"""
function reinterpret_time(time_string)

    split_time = split(time_string, ' ')

    times = [uparse(split_time[i]) for i = 2:length(split_time)]
    sum_time = sum(times) |> u"s" |> ustrip

    return Float64(sum_time)
end


"""
    parse_timing(filename)

Reads the output log file and returns a tuple of Arrays of `(step_number, timing, active_p)`.
"""
function parse_timing(filename::String, Nspecies::Integer=2, output_file::String=""; verbose::Bool=true)

    # read the balance file
    f = open(filename)
    lines = readlines(f)
    close(f)

    # select all lines that contain step information
    sel = findall(occursin.("Step:", lines))[2:end-1]

    # read step counter 
    steps = [parse(Int64, lines[sel[i]][6:15]) for i = 1:length(sel)]

    # read number of particles per species
    active = Vector{Vector{Float64}}(undef, Nspecies)
    for species = 1:Nspecies
        active[species] = [parse(Float64, lines[sel[i]+18+species][30:37]) for i = 1:length(sel)]
    end

    sum_active = Vector{Float64}(undef, length(sel))
    for i = 1:length(sel)
        sum_active[i] = 0.0
        for j = 1:Nspecies
            sum_active[i] +=  active[j][i]
        end
    end

    # get time of step
    # select all lines that contain step information
    timing = [reinterpret_time(lines[sel[i]+20+Nspecies][19:end]) for i = 1:length(sel)]

    if output_file == ""
        return steps, timing, sum_active
    else
        open(output_file; write=true) do f
            write(f, "# Nstep timing active_p\n")
            writedlm(output_file, hcat(steps, timing, sum_active), ' ')
        end
    end
end

