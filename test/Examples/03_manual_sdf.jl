using stl2sdf
using stl2sdf.TerminalUtils
using stl2sdf.DataImport
using stl2sdf.MeshGrid
using stl2sdf.SignedDistances
using stl2sdf.SdfSmoothing
using stl2sdf.DataExport

# ============================================================================
# Manual SDF Conversion Pipeline
# This example shows the step-by-step process without using stl_to_sdf()
# ============================================================================

# Configuration
taskName = "artifacts"
remove_artifacts = true   # Remove small disconnected interior components
artifact_ratio = 0.05     # Keep components ≥5% of largest component
smoothing_method = :interpolation  # :interpolation (exact values) or :approximation (smoother)
grid_refinement = 1       # 1 = same resolution, 2 = double resolution

# 1. Import STL file
# Loads triangular mesh from STL (ASCII or binary format)
print_info("Importing STL file: $(taskName).stl")
(X, IEN) = import_stl("data/$(taskName).stl")

# 2. Create triangular mesh structure
# Converts raw vertex and connectivity data into TriangularMesh type
print_info("Creating triangular mesh")
TriMesh = TriangularMesh(X, IEN)

# 3. Setup SDF computational grid
# Calculate grid dimensions based on mesh bounding box and desired cell size
print_info("Setting up SDF grid")
sdf_grid = interactive_sdf_grid_setup(TriMesh)
points = generateGridPoints(sdf_grid)  # Generate actual 3D grid points

# 4. Compute unsigned distance field
# Calculate shortest distance from each grid point to mesh surface
print_info("Computing unsigned distances")
(dists, xp) = evalDistancesOnTriMesh(TriMesh, sdf_grid, points)

# 5. Determine inside/outside signs using raycasting
# Positive = inside mesh, negative = outside mesh
print_info("Computing signs via raycasting")
(signs, confidences) = raycast_sign_detection(TriMesh, sdf_grid, points)

# Check for low confidence points (may indicate mesh issues)
if any(c -> c < 0.6, confidences)
    low_conf_count = count(c -> c < 0.6, confidences)
    print_warning("$(low_conf_count) points have low confidence (<0.6)")
end

# 6. Combine distances and signs to create signed distance field
print_info("Combining distances and signs to create SDF")
sdf_dists = dists .* signs

# 7. Remove small disconnected components (artifacts)
# Optional step to clean up numerical noise and small floating geometry
if remove_artifacts
    print_info("Removing SDF artifacts")
    nodes_flipped =
        remove_sdf_artifacts!(sdf_dists, sdf_grid, min_component_ratio = artifact_ratio)
    print_success("Artifact removal: $nodes_flipped nodes modified")
end

# 8. Apply RBF smoothing for high-quality SDF
# Interpolation preserves exact values at grid points
# Approximation creates smoother but less accurate results
print_info("Applying RBF smoothing")
is_interpolation = (smoothing_method == :interpolation)
(fine_sdf, fine_grid) =
    RBFs_smoothing(sdf_dists, sdf_grid, is_interpolation, grid_refinement)

# 9. Export results to VTI format for visualization in ParaView
print_info("Exporting results")
exportSdfToVTI("$(taskName)_raw_sdf.vti", sdf_grid, sdf_dists, "distance")
exportSdfToVTI("$(taskName)_smooth_sdf.vti", fine_grid, fine_sdf, "distance")

print_success("SDF conversion completed!")
