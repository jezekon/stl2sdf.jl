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

# Compute local coords:
include("sdfOnTriangularMesh.jl")
include("SignDetection.jl")

end
