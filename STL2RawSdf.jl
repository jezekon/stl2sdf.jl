using stl2sdf
using Dates

taskName = "3D_2x1x1_4Legs_16tol_r2.0"

options = Options(
    cell_size = 0.05,
    smoothing_method = nothing,
    remove_artifacts = false,
    export_raw_sdf = true,
)

# Measure time
start_time = now()
elapsed = @elapsed stl_to_sdf("$(taskName)-STL.stl", options = options)

# Write summary
open("$(taskName)_summary.txt", "w") do f
    println(f, "==================================================")
    println(f, "STL TO SDF CONVERSION SUMMARY")
    println(f, "==================================================")
    println(f, "Task name:           $taskName")
    println(f, "Cell size:           $(options.cell_size)")
    println(f, "Total time:          $(round(elapsed, digits=2)) s")
    println(f, "Generated:           $(Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))")
    println(f, "==================================================")
end
