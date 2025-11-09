using stl2sdf

# Complex SDF conversion with all options specified
# Demonstrates artifact removal and custom grid settings
options = Options(
    smoothing_method = :interpolation,  # Preserves exact SDF values at grid points (vs :approximation for smoother results)
    grid_refinement = 1,                # Output grid resolution: 1=same as input, 2=double resolution
    cell_size = 0.5,                    # Grid spacing in spatial units (nothing=interactive prompt)
    remove_artifacts = true,            # Remove small disconnected interior components
    artifact_ratio = 0.05,              # Keep components ≥5% of largest component size
)

stl_to_sdf("data/artifacts.stl", options = options)
