using Test
using LinearAlgebra
using stl2sdf
using stl2sdf.MeshGrid
# using stl2sdf.SignedDistances
# using stl2sdf.SdfSmoothing

@testset "stl2sdf.jl" begin
    # Test configuration flags
    RUN_lin_beam = true

    if RUN_lin_beam
      taskName = "beam-approx"

      # # Data from Matlab:
      data = read("../data/$(taskName).stl")
      # data = matread("test/" * taskName * ".mat")
      
      (X, IEN) = GenerateMesh.MeshInformations(data)

      # # input data propertis (mesh, density)
      mesh = Rho2sdf.Mesh(X, IEN)
    
    end
end
