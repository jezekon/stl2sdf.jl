module MeshGrid


export MeshInformations, Mesh, TriangularMesh, Grid, LinkedList, getMesh_AABB, generateGridPoints, extractSurfaceTriangularMesh, calculateMiniAABB_grid, generateConnectivityArray, find_triangle_position, interactive_sdf_grid_setup, noninteractive_sdf_grid_setup, NodePosition3D

using Statistics
using LinearAlgebra
using StaticArrays
using Printf
using Base.Threads
using stl2sdf.TerminalUtils

include("MeshVolume.jl")
include("MeshInformations.jl")
include("MeshUtils.jl")
include("SurfaceTriangularMesh.jl")
include("Grid.jl")
include("Grid_setup.jl")


end

