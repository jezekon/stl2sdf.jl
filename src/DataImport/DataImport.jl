module DataImport
export import_tetgen_mesh, import_stl, run_tetgen
 
using LinearAlgebra
using stl2sdf.TerminalUtils

# Import data from STL:
include("ImportStlData.jl")

# Processing and importing Tetgen data:
include("TetgenProcessor.jl")

# Import data from Tetgen (.1):
include("ImportTetgenData.jl")

end
