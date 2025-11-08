using Test
using LinearAlgebra
using stl2sdf
using stl2sdf.TerminalUtils
using stl2sdf.DataImport
using stl2sdf.MeshGrid
using stl2sdf.SignedDistances
using stl2sdf.SdfSmoothing
using stl2sdf.DataExport
using BenchmarkTools

@testset "stl2sdf.jl" begin
    # Test configuration flags
    Import_beam = true
    RUN_main_opts = true
    RUN_beam = true

    if Import_beam
        taskName = "Beam"

        @time (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors
        println(typeof(X))
        println(size(X))
        println(X[1])
        println(typeof(IEN))
        println(size(IEN))
        println(IEN[1])
    end

    if RUN_main_opts
        # Specify all options
        options = Options(
            smoothing_method = :interpolation,
            grid_refinement = 2,
            cell_size = 0.8,
            remove_artifacts = true,
            artifact_ratio = 0.01,
        )
        result = stl_to_sdf("../data/Beam.stl", options = options)
    end

    if RUN_beam
        taskName = "Beam"
        N = 80  # Number of cells along the longest side

        # 1. Import STL file
        print_info("Importing STL file: $taskName")
        (X, IEN) = import_stl("../data/$(taskName).stl")

        # 2. Create triangular mesh (always needed)
        print_info("Creating triangular mesh")
        TriMesh = TriangularMesh(X, IEN)

        # 3. Create SDF grid
        print_info("Setting up SDF grid")
        X_min, X_max = getMesh_AABB(TriMesh.X)
        sdf_grid = Grid(X_min, X_max, N, 3)
        points = generateGridPoints(sdf_grid)

        # 4. Compute SDF
        print_info("Computing unsigned distances")
        (dists, xp) = evalDistancesOnTriMesh(TriMesh, sdf_grid, points)

        print_info("Computing signs")
        (signs, confidences) = raycast_sign_detection(TriMesh, sdf_grid, points)

        # Low confidence points
        if any(c -> c < 0.6, confidences)
            low_conf_count = count(c -> c < 0.6, confidences)
            print_warning("$(low_conf_count) points have low confidence (<0.6)")
        end

        print_info("Combining distances and signs to create SDF")
        sdf_dists = dists .* signs

        # 7. Apply RBF smoothing
        print_info("Applying RBF smoothing")
        (fine_sdf, fine_grid) = RBFs_smoothing(sdf_dists, sdf_grid, true, 1)

        exportSdfToVTI("$(taskName)_sdf.vti", sdf_grid, sdf_dists, "distance")
        exportSdfToVTI("fine-$(taskName)_sdf.vti", fine_grid, fine_sdf, "distance")
    end
end
