using Test
using EntityTools # Replace with your package's actual module name if different
using Makie

@testset "EntityTools.jl" begin

    @testset "Sampling Methods" begin
        @testset "Monopole Sampling" begin
            nth = 30
            mono = monopole_sampling(nth=nth)
            
            @test length(mono) == nth
            @test all(mono .> 0.0)
            @test all(mono .< π)
            @test issorted(mono)
        end

        @testset "Dipole Sampling" begin
            nth = 30
            pole = 1/16
            dip = dipole_sampling(nth=nth, pole=pole)
            
            # Math logic from the Python translation
            expected_poles = floor(Int, nth * pole)
            expected_equator = div(nth - 2 * expected_poles, 2)
            expected_length = 2 * expected_poles + expected_equator
            
            @test length(dip) == expected_length
            @test all(dip .> 0.0)
            @test all(dip .< π)
            @test issorted(dip)
        end
    end

    @testset "Fieldline Integration" begin
        # Create a simple mock grid and uniform vector field
        r_coords = collect(1.0:0.5:5.0)
        th_coords = collect(0.0:0.2:π)
        
        fr_data = ones(length(th_coords), length(r_coords))
        fth_data = zeros(length(th_coords), length(r_coords))
        
        start_pts = [[2.0, π/4], [3.0, π/2]]
        
        lines = compute_fieldlines(r_coords, th_coords, fr_data, fth_data, start_pts, maxsteps=10)
        
        @test length(lines) == 2
        @test length(lines[1]) > 1
        @test length(lines[1][1]) == 2 # Output should be [x, y] arrays
        @test typeof(lines) <: AbstractVector
    end

    @testset "Plotting Recipes (Makie)" begin
        # Test that plotting functions execute without errors
        fig = Figure()
        ax = Axis(fig[1, 1])
        
        r_coords = collect(1.0:1.0:5.0)
        th_coords = collect(0.0:0.5:π)
        data = rand(length(th_coords), length(r_coords))
        
        fr_data = rand(length(th_coords), length(r_coords))
        fth_data = rand(length(th_coords), length(r_coords))
        
        @test_nowarn polar_pcolor!(ax, r_coords, th_coords, data)
        @test_nowarn polar_contour!(ax, r_coords, th_coords, data)
        
        # Test fieldlines plot with the dictionary template
        sample_dict = Dict(:template => "monopole", :nth => 5, :radius => 2.0)
        @test_nowarn plot_fieldlines!(ax, r_coords, th_coords, fr_data, fth_data; 
                                      sample_template=sample_dict, maxsteps=5)
    end

end