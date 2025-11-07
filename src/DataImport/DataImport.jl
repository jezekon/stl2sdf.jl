module DataImport

export import_stl

using FileIO  # FileIO provides the load function that uses MeshIO internally
using LinearAlgebra
using MeshIO
using GeometryBasics


using stl2sdf.TerminalUtils

# Import data from STL:
include("ImportStlData.jl")

end
