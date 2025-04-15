"""
   SelectProjectedNodes(mesh, grid, xp, points)

Extracts nodes that have valid projections on the mesh surface.

This function filters grid points based on their projection status, keeping only 
those with non-zero projections. It also computes statistics on projection distances.

# Arguments
- `mesh::Mesh`: Reference mesh data structure
- `grid::Grid`: Grid containing the point distribution
- `xp::Matrix{Float64}`: Matrix of projected points
- `points::Matrix{Float64}`: Original grid points 

# Returns
- `X`: Vector of points with valid projections
- `Xp`: Vector of corresponding projection points
- `mean_PD`: Mean projection distance
- `max_PD`: Maximum projection distance
"""
function SelectProjectedNodes(
   mesh::Mesh,
   grid::Grid,
   xp::Matrix{Float64},
   points::Matrix{Float64})
   ngp = grid.ngp # number of nodes in grid
   nsd = mesh.nsd # number of spacial dimensions

   # Preallocate arrays with maximum possible size
   max_size = ngp * 2  # Adjust this based on your knowledge of the data
   X = [zeros(Float64, nsd) for _ in 1:max_size]
   Xp = [zeros(Float64, nsd) for _ in 1:max_size]

   count = 0
   for i = 1:ngp
       if sum(abs.(xp[:, i])) > 1.0e-10
           count += 1
           X[count] = points[:, i]
           Xp[count] = xp[:, i]
       end
   end

   # Trim the unused preallocated space
   X = resize!(X, count)
   Xp = resize!(Xp, count)

   # Mean and max projected distance:
   if count > 0
       mean_PD = mean(norm.(X-Xp))
       max_PD = maximum(norm.(X-Xp))
   else 
       mean_PD = 0.
       max_PD = 0.
   end

   return X, Xp, mean_PD, max_PD
end
