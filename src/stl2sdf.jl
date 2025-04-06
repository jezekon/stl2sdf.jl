module stl2sdf

export MeshInformations

using LinearAlgebra
using Statistics
using DelimitedFiles
using BenchmarkTools

# Tools for monitoring the computation process:
include("TerminalUtils/TerminalUtils.jl")
using .TerminalUtils

# Import data from stl and Tetgen results:
include("DataImport/DataImport.jl")
using .DataImport

include("MeshGrid/MeshGrid.jl")
using .MeshGrid

include("DataExport/DataExport.jl")
using .DataExport

include("SignedDistances/SignedDistances.jl")
using .SignedDistances

include("SdfSmoothing/SdfSmoothing.jl")
using .SdfSmoothing

# include("ImplicitDomainMeshing/ImplicitDomainMeshing.jl")
# using .ImplicitDomainMeshing

export import_tetgen_mesh, import_stl

end # module stl2sdf
