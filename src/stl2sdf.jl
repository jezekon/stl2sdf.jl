module stl2sdf

export MeshInformations

using LinearAlgebra
using Statistics
using JLD2

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
include("stl_to_sdf.jl")

export stl_to_sdf, Options

end # module stl2sdf
