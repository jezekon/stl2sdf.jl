# Define a struct to hold optional parameters
struct SDFOptions
    smoothing_method::Symbol
    grid_refinement::Int
    grid_step::Union{Float64,Nothing}

    # Constructor with default values
    function SDFOptions(;
        smoothing_method::Symbol = :interpolation,
        grid_refinement::Int = 1,
        grid_step::Union{Float64,Nothing} = nothing,
    )
        return new(smoothing_method, grid_refinement, grid_step)
    end
end

function stl_to_sdf(stl_filename::String; options::SDFOptions = SDFOptions())
    # Check inputs
    @assert options.smoothing_method in [:interpolation, :approximation] "Smoothing method must be :interpolation or :approximation"
    @assert options.grid_refinement in [1, 2] "Grid refinement must be 1 or 2"

    # Extract base name without extension
    base_name = splitext(basename(stl_filename))[1]

    # 1. Import STL file
    print_info("Importing STL file: $stl_filename")
    (X, IEN) = import_stl(stl_filename)

    # 2. Try Tetgen, fallback to raycast if fails
    TetMesh = nothing
    # 3. Create triangular mesh (always needed)
    print_info("Creating triangular mesh")
    TriMesh = TriangularMesh(X, IEN)

    # 4. Create SDF grid
    print_info("Setting up SDF grid")
    if options.grid_step === nothing
        sdf_grid = interactive_sdf_grid_setup(TriMesh)
    else
        sdf_grid = noninteractive_sdf_grid_setup(TriMesh, options.grid_step)
    end
    points = generateGridPoints(sdf_grid)

    # 5. Compute SDF
    print_info("Computing unsigned distances")
    (dists, xp) = evalDistancesOnTriMesh(TriMesh, sdf_grid, points)

    print_info("Computing signs")
    (signs, confidences) = raycast_sign_detection(TriMesh, sdf_grid, points)

    # Optional: log low confidence points
    if any(c -> c < 0.6, confidences)
        low_conf_count = count(c -> c < 0.6, confidences)
        print_warning("$(low_conf_count) points have low confidence (<0.6)")
    end

    print_info("Combining distances and signs to create SDF")
    sdf_dists = dists .* signs

    # 6. Apply RBF smoothing
    if tetgen_success && TetMesh !== nothing
        print_info("Applying RBF smoothing")
        is_interpolation = (options.smoothing_method === :interpolation)
        (fine_sdf, fine_grid) = RBFs_smoothing(
            TetMesh,
            sdf_dists,
            sdf_grid,
            is_interpolation,
            options.grid_refinement,
            base_name,
        )
    end

    # 7. Save data to JLD2 file
    # sdf_output_file = "$(base_name)_sdf.jld2"
    # grid_output_file = "$(base_name)_grid.jld2"
    # print_info("Saving SDF data to $sdf_output_file and $grid_output_file")
    # @save sdf_output_file fine_sdf
    # @save grid_output_file fine_grid

    # Export the results as VTI format for visualization
    print_info("Exporting SDF to VTI file")
    exportSdfToVTI("$(base_name)_sdf.vti", sdf_grid, sdf_dists, "distance")

    # Return results
    # return (sdf_dists, sdf_grid, fine_sdf, fine_grid)
end
