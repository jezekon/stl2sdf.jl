# src/stl_to_sdf.jl

struct Options
    smoothing_method::Union{Symbol,Nothing}
    grid_refinement::Int
    cell_size::Union{Float64,Nothing}
    remove_artifacts::Bool
    artifact_ratio::Float64
    export_raw_sdf::Bool  # New: export raw SDF to VTI + JLD2

    function Options(;
        smoothing_method::Union{Symbol,Nothing} = nothing,
        grid_refinement::Int = 1,
        cell_size::Union{Float64,Nothing} = nothing,
        remove_artifacts::Bool = false,
        artifact_ratio::Float64 = 0.01,
        export_raw_sdf::Bool = true,
    )
        return new(
            smoothing_method,
            grid_refinement,
            cell_size,
            remove_artifacts,
            artifact_ratio,
            export_raw_sdf,
        )
    end
end

function stl_to_sdf(stl_filename::String; options::Options = Options())
    # Validate inputs
    if options.smoothing_method !== nothing
        @assert options.smoothing_method in [:interpolation, :approximation] "smoothing_method must be :interpolation, :approximation, or nothing"
    end
    @assert options.grid_refinement in [1, 2] "grid_refinement must be 1 or 2"
    @assert options.artifact_ratio > 0.0 && options.artifact_ratio < 1.0 "artifact_ratio must be in range (0, 1)"

    # Extract base name without extension
    base_name = splitext(basename(stl_filename))[1]

    # 1. Import STL file
    print_info("Importing STL file: $stl_filename")
    (X, IEN) = import_stl(stl_filename)

    # 2. Create triangular mesh
    print_info("Creating triangular mesh")
    TriMesh = TriangularMesh(X, IEN)

    # 3. Setup SDF grid
    print_info("Setting up SDF grid")
    if options.cell_size === nothing
        sdf_grid = interactive_sdf_grid_setup(TriMesh)
    else
        sdf_grid = noninteractive_sdf_grid_setup(TriMesh, options.cell_size)
    end
    points = generateGridPoints(sdf_grid)

    # 4. Compute unsigned distances
    print_info("Computing unsigned distances")
    (dists, xp) = evalDistancesOnTriMesh(TriMesh, sdf_grid, points)

    # 5. Determine signs via raycasting
    print_info("Computing signs")
    (signs, confidences) = raycast_sign_detection(TriMesh, sdf_grid, points)

    if any(c -> c < 0.6, confidences)
        low_conf_count = count(c -> c < 0.6, confidences)
        print_warning("$(low_conf_count) points have low confidence (<0.6)")
    end

    # 6. Combine to create signed distance field
    print_info("Combining distances and signs to create SDF")
    sdf_dists = dists .* signs

    # 7. Remove artifacts if enabled
    if options.remove_artifacts
        print_info("Removing SDF artifacts")
        nodes_flipped = remove_sdf_artifacts!(
            sdf_dists,
            sdf_grid,
            min_component_ratio = options.artifact_ratio,
        )
        print_success("Artifact removal completed: $nodes_flipped nodes modified")
    end

    # 8. Export raw SDF if requested
    if options.export_raw_sdf
        B = round(sdf_grid.cell_size, digits = 4)

        # Export to VTI
        exportSdfToVTI(
            "$(base_name)_SDF_CellSize-$(B).vti",
            sdf_grid,
            sdf_dists,
            "distance",
        )

        # Export to JLD2
        println("Saving raw SDF results to JLD2 files...")
        sdf_dists_processed = SdfSmoothing.process_vector(sdf_dists)
        raw_sdf_dims = sdf_grid.N .+ 1
        raw_sdf = reshape(sdf_dists_processed, Tuple(raw_sdf_dims))
        raw_grid = SdfSmoothing.create_grid(sdf_grid.N, sdf_grid)

        @save "Z_$(base_name)_RawSDF_B-$(B).jld2" raw_sdf
        @save "Z_$(base_name)_RawGrid_B-$(B).jld2" raw_grid

        print_success("Raw SDF exported to VTI and JLD2")
    end

    # 9. Apply RBF smoothing (only if smoothing_method is specified)
    if options.smoothing_method !== nothing
        print_info("Applying RBF smoothing")
        is_interpolation = (options.smoothing_method === :interpolation)
        (fine_sdf, fine_grid) =
            RBFs_smoothing(sdf_dists, sdf_grid, is_interpolation, options.grid_refinement)

        exportSdfToVTI("$(base_name)_fine_sdf.vti", fine_grid, fine_sdf, "distance")
    end

    return ()
end
