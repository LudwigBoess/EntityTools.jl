"""
    phase_map(x_q, y_q, x_lim, y_lim; 
              xbins::Int=400, ybins::Int=50,
              show_progress::Bool=true)

Get a 2D histogram of `x_q` and `y_q` in the limits `x_lim`, `y_lim` over a number if bins `xbins` and `ybins`.
"""
function phase_map(x_q, y_q, x_lim, y_lim; 
                   xbins::Int=400, ybins::Int=50,
                   show_progress::Bool=true)

    # get bin spacing
    dx = (x_lim[2] - x_lim[1] ) / xbins
    dy = (y_lim[2] - y_lim[1] ) / ybins

    # allocate ybins x xbins matrix filled with zeros
    phase_map_count = zeros(Int64, xbins, ybins)

    # optional progress meter
    if show_progress
        P = Progress(size(x_q,1))
        idx_p = 0
    end

    @inbounds for i = 1:size(x_q,1)

        x_bin = 1 + floor( Int64, (x_q[i] - x_lim[1])/dx )
        y_bin = 1 + floor( Int64, (y_q[i] - y_lim[1])/dy )

        if (1 <= x_bin <= xbins) && (1 <= y_bin <= ybins)
            phase_map_count[x_bin, y_bin] += 1   
        end

        # update progress meter
        if show_progress
            idx_p += 1
            ProgressMeter.update!(P, idx_p)
        end
    end

    # return 2D histogram
    return LinRange(x_lim[1], x_lim[2], xbins), # X
           LinRange(y_lim[1], y_lim[2], ybins), # Y
           phase_map_count
end