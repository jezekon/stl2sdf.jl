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
    Import_beam = false
    Import_bunny = false
    RUN_lin_beam = false
    RUN_main_some_opts = false
    RUN_main_opts = false
    RUN_beam_Mecas = true

    #NOTE:
    if Import_beam
        taskName = "beam-approx"

        @time (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors
        println(typeof(X))
        println(size(X))
        println(X[1])
        println(typeof(IEN))
        println(size(IEN))
        println(IEN[1])
    end

    #NOTE:
    if Import_bunny
        # taskName = "StanfordBunny_small"
        taskName = "StanfordBunny_large"

        @time (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors
    end

    #NOTE:
    if RUN_lin_beam
        taskName = "beam-approx"
        # taskName = "StanfordBunny_small"

        N = 120  # Number of cells along the longest side

        # Data from stl:
        (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors
        # (X, IEN) = import_stl("data/$(taskName).stl")

        # Data from Tetgen
        run_tetgen(taskName, "../data")  # Run in specified directory
        (X_tet, IEN_tet) = import_tetgen_mesh("../data/$(taskName).1") # -> Vector of vectors
        # (X_tet, IEN_tet) = import_tetgen_mesh("data/$(taskName).1")

        TriMesh = TriangularMesh(X, IEN)
        TetMesh = Mesh(X_tet, IEN_tet)

        ## Grid:
        X_min, X_max = getMesh_AABB(TriMesh.X)
        sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
        points = generateGridPoints(sdf_grid) # uzly pravidelné mřížky

        ## SFD from triangular mesh:
        (dists, xp) = evalDistancesOnTriMesh(TriMesh, sdf_grid, points) # Vector{Float64}
        (signs, _) = raycast_sign_detection(TriMesh, sdf_grid, points)
        sdf_dists = dists .* signs

        exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

        # RBF smoothing:
        is_interp = false
        (fine_sdf, fine_grid) =
            RBFs_smoothing(TetMesh, sdf_dists, sdf_grid, is_interp, 2, taskName) # interpolation == true, aproximation == false, smooth
        export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, 2, is_interp)
    end

    if RUN_main_some_opts
        options = SDFOptions(grid_step = 0.8)  # Use defaults for other parameters
        result = stl_to_sdf("beam-approx.stl", options = options)
    end

    if RUN_main_opts
        # Specify all options
        options = SDFOptions(
            smoothing_method = :interpolation,
            grid_refinement = 2,
            grid_step = 0.8,
        )
        result = stl_to_sdf("beam-approx.stl", options = options)
    end

    if RUN_beam_Mecas
        taskName = "beam_Mecas"
        N = 80  # Number of cells along the longest side

        # 1. Import STL file
        print_info("Importing STL file: $taskName")
        (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors

        # 3. Create triangular mesh (always needed)
        print_info("Creating triangular mesh")
        TriMesh = TriangularMesh(X, IEN)

        # 4. Create SDF grid
        print_info("Setting up SDF grid")
        X_min, X_max = getMesh_AABB(TriMesh.X)
        sdf_grid = Grid(X_min, X_max, N, 3) # cartesian grid
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
        print_info("Applying RBF smoothing")
        (fine_sdf, fine_grid) = RBFs_smoothing(
            sdf_dists,
            sdf_grid,
            is_interpolation = :interpolation,
            grid_refinement = 1.0,
        )

        exportSdfToVTI("$(taskName)_sdf.vti", sdf_grid, sdf_dists, "distance")
        exportSdfToVTI("fine-$(base_name)_sdf.vti", fine_grid, fine_sdf, "distance")
    end
end
