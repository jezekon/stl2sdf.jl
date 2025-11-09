# src/stl_to_sdf.jl

"""
    Options

Configuration for STL to SDF conversion process.

# Fields
- `smoothing_method::Union{Symbol,Nothing}` - RBF smoothing method: :interpolation, :approximation, or nothing (no smoothing)
- `grid_refinement::Int` - Grid refinement factor for smoothed SDF (1 = same, 2 = double resolution)
- `cell_size::Union{Float64,Nothing}` - Grid cell size in spatial units (nothing triggers interactive setup)
- `remove_artifacts::Bool` - Remove small disconnected components from SDF field
- `artifact_ratio::Float64` - Minimum component size as ratio of largest component (0.01 = 1%)

# Example
```julia
opts = Options(cell_size=0.5, remove_artifacts=true, smoothing_method=nothing)  # No smoothing
opts = Options(cell_size=0.5, smoothing_method=:interpolation, grid_refinement=2)  # With smoothing
```
"""
struct Options
    smoothing_method::Union{Symbol,Nothing}
    grid_refinement::Int
    cell_size::Union{Float64,Nothing}
    remove_artifacts::Bool
    artifact_ratio::Float64

    function Options(;
        smoothing_method::Union{Symbol,Nothing} = nothing,
        grid_refinement::Int = 1,
        cell_size::Union{Float64,Nothing} = nothing,
        remove_artifacts::Bool = false,
        artifact_ratio::Float64 = 0.01,
    )
        return new(
            smoothing_method,
            grid_refinement,
            cell_size,
            remove_artifacts,
            artifact_ratio,
        )
    end
end

"""
    stl_to_sdf(stl_filename::String; options::Options = Options())

Convert STL mesh to signed distance function with optional RBF smoothing.

# Arguments
- `stl_filename::String` - Path to input STL file (ASCII or binary format)
- `options::Options` - Configuration options for SDF generation

# Returns
- Nothing (exports VTI files)

# Workflow
1. Import and process STL mesh
2. Setup computational grid (interactive or automatic)
3. Compute unsigned distances to mesh surface
4. Determine signs via raycasting
5. Optionally remove artifacts from SDF
6. Optionally apply RBF smoothing with optional refinement
7. Export VTI files for visualization

# Example
```julia
opts = Options(cell_size=0.5, remove_artifacts=true, smoothing_method=nothing)
stl_to_sdf("model.stl", options=opts)  # No smoothing, only raw SDF

opts = Options(cell_size=0.5, smoothing_method=:interpolation, grid_refinement=2)
stl_to_sdf("model.stl", options=opts)  # With smoothing
```
"""
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

    # Log low confidence points
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

    # 8. Export raw SDF result
    exportSdfToVTI("$(base_name)_sdf.vti", sdf_grid, sdf_dists, "distance")

    # 9. Apply RBF smoothing (only if smoothing_method is specified)
    if options.smoothing_method !== nothing
        print_info("Applying RBF smoothing")
        is_interpolation = (options.smoothing_method === :interpolation)
        (fine_sdf, fine_grid) =
            RBFs_smoothing(sdf_dists, sdf_grid, is_interpolation, options.grid_refinement)

        # Export smoothed result
        exportSdfToVTI("$(base_name)_fine_sdf.vti", fine_grid, fine_sdf, "distance")
    end

    return ()
end
