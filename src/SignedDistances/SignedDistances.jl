module SignedDistances

export evalDistancesOnTriMesh, raycast_sign_detection

using Base.Threads
using Statistics
using StaticArrays
using LinearAlgebra
using ProgressMeter
using NearestNeighbors
using Random
using ImplicitBVH
using ImplicitBVH: BBox, BSphere

using stl2sdf.TerminalUtils
using stl2sdf.MeshGrid
using stl2sdf.DataExport
using stl2sdf

# Compute the unsigned distance field:
include("dfOnTriangularMesh.jl")

# Raycast sign detection:
include("RaycastingSignDetection.jl")

end
