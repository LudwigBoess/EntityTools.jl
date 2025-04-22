"""
    get_histograms(bins, Q)

Computes the histograms for `M_s` of all shocked particles
"""
function get_histograms(bins, Q)

    # fit histograms
    hist = fit(Histogram, Q, bins)

    return hist.weights
end

"""
    construct_bins_and_centers(bin_min, bin_max, Nbins)

Returns the bins and their center points.
"""
function construct_bins_and_centers(bin_min, bin_max, Nbins)

    # construct bin boundaries
    bins = LinRange(bin_min, bin_max, Nbins + 1)

    # construct bin centers
    bin_centers = Vector{Float64}(undef, Nbins)

    for i = 1:Nbins
        bin_centers[i] = 0.5 * (bins[i] + bins[i+1])
    end

    return bins, bin_centers
end


"""
    spectrum(Q, spec_min, spec_max; Nbins=100)

Constructs a spectrum `dQ/dx` between x values between `spec_min` and `spec_max`.
Uses a number of `Nbins=100` by default.
"""
function spectrum(Q, spec_min, spec_max; Nbins=100)

    bins, bin_centers = construct_bins_and_centers(log10(spec_min), log10(spec_max), Nbins)

    bin_centers = 10.0 .^ bin_centers
    bins = 10.0 .^ bins

    dQ = [bins[i+1] - bins[i] for i = 1:length(bins)-1]

    Q_hist = get_histograms(bins, Q)

    return bin_centers, Q_hist ./ dQ
end