module DataImport
export import_tetgen_mesh, import_stl
 
using LinearAlgebra

# Import data from STL:
include("ImportStlData.jl")

# Import data from Tetgen (.1):
include("ImportTetgenData.jl")

end
