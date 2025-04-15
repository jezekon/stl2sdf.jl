module SignedDistances

export evalDistancesOnTriMesh, SignDetection

using Base.Threads
using Statistics
using StaticArrays
using LinearAlgebra
using ProgressMeter
using BenchmarkTools
using NearestNeighbors

using stl2sdf.TerminalUtils
using stl2sdf.MeshGrid
using stl2sdf.DataExport
using stl2sdf

# Compute the unsigned distance field:
include("dfOnTriangularMesh.jl")

# Compute signs:
include("SignDetection.jl")

end
