module SignedDistances

export evalDistancesOnTriMesh, SignDetection, RaycastingSignDetection

using Base.Threads
using Statistics
using StaticArrays
using LinearAlgebra
using ProgressMeter
using NearestNeighbors
using Random

using stl2sdf.TerminalUtils
using stl2sdf.MeshGrid
using stl2sdf.DataExport
using stl2sdf

# Compute the unsigned distance field:
include("dfOnTriangularMesh.jl")

# Compute signs:
include("SignDetection.jl")

# Raycast fallback:
include("RaycastingSignDetection.jl")

end
