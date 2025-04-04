module SignedDistances

export evalDistancesOnTriangularMesh, Sign_Detection

using Base.Threads
# using Einsum
using Statistics
using StaticArrays
using LinearAlgebra
# using DelimitedFiles
using ProgressMeter
# using NLopt
using BenchmarkTools

using Rho2sdf.TerminalUtils
using Rho2sdf.MeshGrid
using Rho2sdf

# Compute local coords:
include("PseudoNormals.jl")
include("sdfOnTriangularMesh.jl")
include("SignDetection.jl")

end
