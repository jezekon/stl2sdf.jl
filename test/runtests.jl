using Test
using LinearAlgebra
using stl2sdf
using stl2sdf.DataImport
using stl2sdf.MeshGrid
using stl2sdf.SignedDistances
using stl2sdf.SdfSmoothing
using stl2sdf.DataExport

@testset "stl2sdf.jl" begin
    # Test configuration flags
    Import_beam = false
    Import_bunny = false
    RUN_lin_beam = true
    RUN_main_some_opts = false
    RUN_main_opts = false
    
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
      signs = @time SignDetection(TetMesh, sdf_grid, points)
      sdf_dists = dists .* signs

      exportSdfToVTI(taskName * "_SDF.vti", sdf_grid, sdf_dists, "distance")

      # RBF smoothing:
      is_interp = false
      (fine_sdf, fine_grid) = RBFs_smoothing(TetMesh, sdf_dists, sdf_grid, is_interp, 2, taskName) # interpolation == true, aproximation == false, smooth
      export_sdf_results(fine_sdf, fine_grid, sdf_grid, taskName, 2, is_interp)
    end

    if RUN_main_some_opts
      options = SDFOptions(grid_step = 0.8)  # Use defaults for other parameters
      result = stl_to_sdf("beam-approx.stl", options=options)
    end

    if RUN_main_opts
      # Specify all options
      options = SDFOptions(
          smoothing_method = :interpolation,
          grid_refinement = 2,
          grid_step = 0.8
      )
      result = stl_to_sdf("beam-approx.stl", options=options)
    end
end
