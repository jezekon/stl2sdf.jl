module DataExport
   
export SelectProjectedNodes, exportToVTU, exportStructuredPointsToVTK

using Statistics
using LinearAlgebra
using WriteVTK
using stl2sdf
using stl2sdf.MeshGrid

include("DataPostProcess.jl")
include("ExportToVTU.jl")
include("ExportToVTK.jl")

end
