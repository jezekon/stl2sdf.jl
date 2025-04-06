using Test
using LinearAlgebra
using stl2sdf
using stl2sdf.DataImport
using stl2sdf.MeshGrid
using stl2sdf.SignedDistances
# using stl2sdf.SdfSmoothing
using stl2sdf.DataExport

@testset "stl2sdf.jl" begin
    # Test configuration flags
    RUN_lin_beam = true

    if RUN_lin_beam
      taskName = "beam-approx"
      N = 40  # Number of cells along the longest side

      # Data from stl:
      (X, IEN) = import_stl("../data/$(taskName).stl") # -> Vector of vectors
      # (X, IEN) = import_stl("data/$(taskName).stl")
      
      # Data from Tetgen
      run_tetgen("beam-approx", "../data")  # Run in specified directory
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

      exportStructuredPointsToVTK(taskName * "_SDF.vtk", sdf_grid, sdf_dists, "distance")

    end
end
