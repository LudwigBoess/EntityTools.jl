"""
    dipole_sampling(; nth::Int=30, pole::Float64=1/16)

Returns an array of angles sampled from a dipole distribution.
"""
function dipole_sampling(; nth::Int=30, pole::Float64=1/16)
    nth_poles = floor(Int, nth * pole)
    nth_equator = div(nth - 2 * nth_poles, 2)
    
    # Generate ranges mimicking np.linspace behavior (omitting endpoints where appropriate)
    p1 = range(0, stop=π * pole, length=nth_poles + 1)[2:end]
    p2 = range(π * pole, stop=π / 2, length=nth_equator + 2)[2:end-1]
    p3 = range(π * (1 - pole), stop=π, length=nth_poles + 1)[1:end-1]
    
    return vcat(p1, p2, p3)
end

"""
    monopole_sampling(; nth::Int=30)

Returns an array of angles sampled from a monopole distribution.
"""
function monopole_sampling(; nth::Int=30)
    return range(0, stop=π, length=nth + 2)[2:end-1]
end

"""
    compute_fieldlines(r_coords, th_coords, fr_data, fth_data, start_points; kwargs...)

Computes field lines of a vector field defined by radial (`fr`) and azimuthal (`fth`) components.
"""
function compute_fieldlines(r_coords::AbstractVector, th_coords::AbstractVector, 
                            fr_data::AbstractMatrix, fth_data::AbstractMatrix, 
                            start_points::AbstractVector;
                            direction::String="both", 
                            stop_when::Function=(xy, rth)->false, 
                            ds::Float64=0.1, 
                            maxsteps::Int=1000)
    
    # Calculate vector field components in Cartesian space for integration
    # Note: Assuming data matrices are ordered (th, r) to match original code structure
    fxs = zeros(length(th_coords), length(r_coords))
    fys = zeros(length(th_coords), length(r_coords))
    
    for (j, r) in enumerate(r_coords)
        for (i, th) in enumerate(th_coords)
            fxs[i, j] = fr_data[i, j] * sin(th) + fth_data[i, j] * cos(th)
            fys[i, j] = fr_data[i, j] * cos(th) - fth_data[i, j] * sin(th)
        end
    end

    # Nearest neighbor interpolation with 0 fill value
    interp_fx = extrapolate(interpolate((th_coords, r_coords), fxs, Gridded(Constant())), 0.0)
    interp_fy = extrapolate(interpolate((th_coords, r_coords), fys, Gridded(Constant())), 0.0)
    
    rmin, rmax = minimum(r_coords), maximum(r_coords)

    # Internal integration step
    function _integrate(r_th_start, delta)
        r0, th0 = r_th_start
        xy = [r0 * sin(th0), r0 * cos(th0)]
        rth = [r0, th0]
        fieldline = [xy]
        
        for _ in 1:maxsteps
            x, y = xy[1], xy[2]
            r = sqrt(x^2 + y^2)
            th = atan(-y, x) + π / 2
            rth = [r, th]
            
            vx = interp_fx(th, r)
            vy = interp_fy(th, r)
            vmag = sqrt(vx^2 + vy^2)
            
            if vmag == 0 || isnan(vmag)
                break
            end
            
            xy = xy .+ delta .* [vx, vy] ./ vmag
            
            stop_condition = stop_when(xy, rth) || (rth[1] < rmin) || (rth[1] > rmax) || 
                             (rth[2] < 0) || (rth[2] > π)
                             
            if stop_condition || any(isnan.(xy)) || any(isinf.(xy))
                break
            else
                push!(fieldline, xy)
            end
        end
        return fieldline
    end

    lines = []
    for pt in start_points
        if direction == "forward"
            push!(lines, _integrate(pt, ds))
        elseif direction == "backward"
            push!(lines, _integrate(pt, -ds))
        else
            f1 = _integrate(pt, ds)
            f2 = _integrate(pt, -ds)
            # Combine lines, avoiding duplication of the starting point
            push!(lines, vcat(reverse(f2)[1:end-1], f1))
        end
    end
    
    return lines
end


"""
    polar_fieldlines!(ax, r_coords, th_coords, fr_data, fth_data; kwargs...)

Plots field lines on an existing Makie axis.
"""
function polar_fieldlines!(ax::Axis, r_coords, th_coords, fr_data, fth_data;
                          start_points=nothing, 
                          sample_template=nothing, 
                          invert_x::Bool=false, 
                          invert_y::Bool=false,
                          direction::String="both",
                          ds::Float64=0.1,
                          maxsteps::Int=1000,
                          line_kwargs...)
    
    if isnothing(start_points) && isnothing(sample_template)
        throw(ArgumentError("Either start_points or sample_template must be specified"))
    elseif isnothing(start_points)
        radius = get(sample_template, :radius, 1.5)
        template = get(sample_template, :template, "dipole")
        
        if template == "dipole"
            start_points = [[radius, th] for th in dipole_sampling(; sample_template...)]
        elseif template == "monopole"
            start_points = [[radius, th] for th in monopole_sampling(; sample_template...)]
        else
            throw(ArgumentError("Unknown sampling template: $template"))
        end
    end

    lines_data = compute_fieldlines(r_coords, th_coords, fr_data, fth_data, start_points; 
                                    direction=direction, ds=ds, maxsteps=maxsteps)
    
    for fl in lines_data
        fl_matrix = reduce(hcat, fl)' # Convert vector of vectors to matrix
        
        x_vals = invert_x ? -fl_matrix[:, 1] : fl_matrix[:, 1]
        y_vals = invert_y ? -fl_matrix[:, 2] : fl_matrix[:, 2]
        
        lines!(ax, x_vals, y_vals; line_kwargs...)
    end
end

"""
    polar_pcolor!(ax, r_coords, th_coords, data_values; kwargs...)

Plots a pseudocolor plot of 2D polar data onto a rectilinear projection.
Creates a triangular mesh to cleanly handle non-uniform cells.
"""
function polar_pcolor!(ax::Axis, r_coords::AbstractVector, th_coords::AbstractVector, 
                       data_values::AbstractMatrix;
                       invert_x::Bool=false, 
                       invert_y::Bool=false,
                       colormap=:viridis,
                       mesh_kwargs...)
    
    nr = length(r_coords)
    nth = length(th_coords)
    
    vertices = Point2f[]
    faces = GLTriangleFace[]
    colors = Float64[]
    
    # Map polar to rectilinear coordinates
    for (i, th) in enumerate(th_coords)
        for (j, r) in enumerate(r_coords)
            x = r * sin(th)
            y = r * cos(th)
            x = invert_x ? -x : x
            y = invert_y ? -y : y
            push!(vertices, Point2f(x, y))
            push!(colors, data_values[j, i])
        end
    end
    
    # Generate triangulation for the grid
    for j in 1:(nth-1)
        for i in 1:(nr-1)
            idx1 = (j-1)*nr + i
            idx2 = (j-1)*nr + i + 1
            idx3 = j*nr + i + 1
            idx4 = j*nr + i
            
            push!(faces, GLTriangleFace(idx1, idx2, idx3))
            push!(faces, GLTriangleFace(idx1, idx3, idx4))
        end
    end
    
    ax.aspect = DataAspect()
    # Use standard Makie mesh rendering for tripped-pseudocolor mapping
    mesh!(ax, vertices, faces, color=colors, colormap=colormap; shading=NoShading, mesh_kwargs...)
end

"""
    polar_contour!(ax, r_coords, th_coords, data_values; kwargs...)

Plots standard contours for 2D polar data on rectilinear axes.
"""
function polar_contour!(ax::Axis, r_coords::AbstractVector, th_coords::AbstractVector, data_values::AbstractMatrix;
                        invert_x::Bool=false, 
                        invert_y::Bool=false,
                        contour_kwargs...)
    
    # Create mapped X and Y grids
    X = zeros(length(th_coords), length(r_coords))
    Y = zeros(length(th_coords), length(r_coords))
    
    for (j, r) in enumerate(r_coords)
        for (i, th) in enumerate(th_coords)
            x = r * sin(th)
            y = r * cos(th)
            X[i, j] = invert_x ? -x : x
            Y[i, j] = invert_y ? -y : y
        end
    end
    
    ax.aspect = DataAspect()
    contour!(ax, X, Y, data_values; contour_kwargs...)
end