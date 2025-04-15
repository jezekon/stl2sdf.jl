module MeshGrid

export MeshInformations, Mesh, TriangularMesh, Grid, LinkedList, getMesh_AABB, generateGridPoints, extractSurfaceTriangularMesh, calculateMiniAABB_grid, generateConnectivityArray, find_triangle_position, interactive_sdf_grid_setup, noninteractive_sdf_grid_setup, NodePosition3D

using Statistics
using LinearAlgebra
using StaticArrays
using Printf
using Base.Threads
using stl2sdf.TerminalUtils

# Functions for computing mesh volume using tetrahedral elements
include("MeshVolume.jl")
# Core mesh data structures and connectivity information
include("MeshInformations.jl")
# Utility functions for mesh manipulation and analysis
include("MeshUtils.jl")
# Triangular surface mesh extraction and operations
include("SurfaceTriangularMesh.jl")
# Grid generation and spatial acceleration structures
include("Grid.jl")
# Interactive and non-interactive grid configuration tools
include("Grid_setup.jl")

end
