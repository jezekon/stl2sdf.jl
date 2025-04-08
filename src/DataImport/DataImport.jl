module DataImport
export import_tetgen_mesh, import_stl, run_tetgen
 
using FileIO  # FileIO provides the load function that uses MeshIO internally
using LinearAlgebra
using MeshIO
using GeometryBasics


using stl2sdf.TerminalUtils

# Import data from STL:
include("ImportStlData.jl")

# Processing and importing Tetgen data:
include("TetgenProcessor.jl")

# Import data from Tetgen (.1):
include("ImportTetgenData.jl")

end
