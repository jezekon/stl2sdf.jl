module SignedDistances

export evalDistancesOnTriMesh,
    raycast_sign_detection, remove_sdf_artifacts!, analyze_sdf_components

using Base.Threads
using Statistics
using StaticArrays
using LinearAlgebra
using ProgressMeter
using NearestNeighbors
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

# SDF artifact removal:
include("SdfArtifactRemoval.jl")

end
