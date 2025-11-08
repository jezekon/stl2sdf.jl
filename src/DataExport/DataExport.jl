module DataExport

export SelectProjectedNodes, exportSdfToVTI, export_sdf_results

using Statistics
using LinearAlgebra
using WriteVTK
using JLD2
using stl2sdf
using stl2sdf.MeshGrid

# Post-processing utilities for filtering and analyzing projected nodes
include("DataPostProcess.jl")

# Exporting SDF grid data to VTI format (image data visualization)
include("ExportToVTI.jl")

# Exporting complete SDF results with metadata
include("ExportSdfResults.jl")

end
