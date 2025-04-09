# Selection of regular grid points that have been projected:
function SelectProjectedNodes(
  mesh::TriangularMesh,
  grid::Grid,
  xp::Matrix{Float64},
  points::Matrix{Float64})

  ngp = grid.ngp # number of nodes in grid
  nsd = mesh.nsd # number of spacial dimensions

  # Assuming ngp is defined somewhere in your code
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

  # If count is 0, indicating no points were added, handle gracefully
  if count == 0
    # println("WARNING: no projected points!")
    return [], [], NaN, NaN
  end

  # Trim the unused preallocated space
  X = resize!(X, count)
  Xp = resize!(Xp, count)

  # Mean and max projected distance:
  mean_PD = mean(norm.(X - Xp))
  max_PD = maximum(norm.(X - Xp))

  return X, Xp, mean_PD, max_PD
end

function barycentricCoordinates(
  x₁::Vector{Float64}, # coordinates of vertex one of the triangle
  x₂::Vector{Float64}, # coordinates of vertex two of the triangle
  x₃::Vector{Float64}, # coordinates of vertex tree of the triangle
  n::Vector{Float64},  # unit normal to the face of the triangle
  x)  # one node of the grid

    A = [
        (x₁[2]*n[3]-x₁[3]*n[2]) (x₂[2]*n[3]-x₂[3]*n[2]) (x₃[2]*n[3]-x₃[3]*n[2])
        (x₁[3]*n[1]-x₁[1]*n[3]) (x₂[3]*n[1]-x₂[1]*n[3]) (x₃[3]*n[1]-x₃[1]*n[3])
        (x₁[1]*n[2]-x₁[2]*n[1]) (x₂[1]*n[2]-x₂[2]*n[1]) (x₃[1]*n[2]-x₃[2]*n[1])
    ]
    b = [
        x[2] * n[3] - x[3] * n[2],
        x[3] * n[1] - x[1] * n[3],
        x[1] * n[2] - x[2] * n[1],
    ]
    
    n_max, i_max = findmax(abs.(n)) ##???
    A[i_max, :] = [1.0 1.0 1.0]
    b[i_max] = 1.0

    return λ = A \ b # barycentric coordinates
end

function evalDistancesOnTriMesh(mesh::TriangularMesh, grid::Grid, points::Matrix{Float64})
    println("Building accelerated spatial structure...")
    # Convert the triangular mesh into a format suitable for KD-tree operations
    # Preprocessing all triangles to create an efficient searchable structure
    triangle_data = prepare_triangle_data(mesh)
    
    println("Creating KD-tree for mesh vertices...")
    # Create a KD-tree from all mesh vertices for nearest neighbor searches
    vtx_kdtree = KDTree(mesh.X)
    
    # Convert the triangle data to a format suitable for range searches
    triangle_centers = Matrix{Float64}(hcat([triangle_data[i].center for i in 1:length(triangle_data)]...))
    
    println("Creating KD-tree for triangle centers...")
    # Create a KD-tree from triangle centers for faster triangle queries
    tri_kdtree = KDTree(triangle_centers)
    
    ngp = grid.ngp
    nsd = mesh.nsd
    nel = mesh.nel
    
    # Use large negative value for initialization of distances
    big = -1.0e10
    dist = big * ones(Float64, ngp)
    xp = zeros(Float64, nsd, ngp)
    
    # Set search radius based on grid cell size 
    # (can be adjusted based on mesh characteristics)
    search_radius = 2.5 * grid.cell_size
    
    println("Computing distances to triangular mesh...")
    
    # Create thread-local storage for parallel processing
    nthreads = Threads.nthreads()
    dist_local = [fill(big, ngp) for _ in 1:nthreads]
    xp_local = [zeros(Float64, nsd, ngp) for _ in 1:nthreads]
    
    # Create a progress bar
    p = Progress(ngp, 1, "Computing distances: ", 30)
    
    # Use atomic counter for progress tracking
    counter = Atomic{Int}(0)
    update_interval = max(1, div(ngp, 100))
    
    # Process points in batches to improve cache efficiency
    batch_size = 128
    num_batches = ceil(Int, ngp / batch_size)
    
    # Parallel processing of points
    Threads.@threads for batch in 1:num_batches
        tid = Threads.threadid()
        start_idx = (batch - 1) * batch_size + 1
        end_idx = min(batch * batch_size, ngp)
        
        for v in start_idx:end_idx
            point = @view points[:, v]
            
            # First, quickly find nearest triangles in a radius
            # This is much faster than checking all triangles or using grid cells
            tri_indices = inrange(tri_kdtree, point, search_radius)
            
            # If no triangles found in range, try with nearest vertex approach
            if isempty(tri_indices)
                # Find nearest vertices to the point
                idxs, dists = knn(vtx_kdtree, point, 3)
                
                # Collect all triangles containing these vertices
                vtx_tri_indices = Set{Int}()
                for vtx_idx in idxs
                    for el in mesh.INE[vtx_idx]
                        push!(vtx_tri_indices, el)
                    end
                end
                tri_indices = collect(vtx_tri_indices)
            end
            
            # Initialize minimum distance for this point
            min_dist = abs(dist_local[tid][v])
            min_proj = zeros(Float64, nsd)
            
            # Process all candidate triangles
            for el in tri_indices
                # Get triangle data
                tri_data = triangle_data[el]
                
                # Quick AABB check (early rejection)
                if !point_in_box_range(point, tri_data.aabb_min, tri_data.aabb_max, min_dist)
                    continue
                end
                
                # Extract triangle vertices
                x₁, x₂, x₃ = tri_data.vertices
                
                # Triangle normal
                n = tri_data.normal
                
                # Compute detailed point-to-triangle distance
                λ = barycentricCoordinates(x₁, x₂, x₃, n, point)
                
                # Initialize projection point and distance check flag
                xₚ = zeros(Float64, nsd)
                isFaceOrEdge = false
                
                # Check if projection is inside triangle face
                if minimum(λ) >= 0.0
                    # Point projects onto triangle face
                    xₚ = λ[1] * x₁ + λ[2] * x₂ + λ[3] * x₃
                    dist_tmp = norm(point - xₚ)
                    
                    # Update if closer than current minimum
                    if dist_tmp < min_dist
                        min_dist = dist_tmp
                        min_proj = copy(xₚ)
                        isFaceOrEdge = true
                    end
                else
                    # Check projections onto triangle edges
                    for j = 1:3
                        edge = tri_data.edges[j]
                        xᵥ = tri_data.vertices[j]
                        L = tri_data.edge_lengths[j]
                        
                        # Project point onto edge line
                        P = dot(point - xᵥ, edge / L)
                        
                        # Check if projection is on the edge segment
                        if P >= 0 && P <= L
                            xₚ = xᵥ + (edge / L) * P
                            dist_tmp = norm(point - xₚ)
                            
                            # Update if closer than current minimum
                            if dist_tmp < min_dist
                                min_dist = dist_tmp
                                min_proj = copy(xₚ)
                                isFaceOrEdge = true
                            end
                        end
                    end
                    
                    # If no face/edge projection found, check vertices
                    if !isFaceOrEdge
                        # Find closest vertex
                        vertex_dists = [norm(point - v) for v in tri_data.vertices]
                        dist_tmp, idx = findmin(vertex_dists)
                        
                        # Update if closer than current minimum
                        if dist_tmp < min_dist
                            min_dist = dist_tmp
                            min_proj = copy(tri_data.vertices[idx])
                        end
                    end
                end
            end
            
            # Update thread-local distance and projection
            if min_dist < abs(dist_local[tid][v])
                dist_local[tid][v] = min_dist
                xp_local[tid][:, v] = min_proj
            end
            
            # Update progress counter
            count = atomic_add!(counter, 1)
            if count % update_interval == 0 && Threads.threadid() == 1
                update!(p, min(count, ngp))
            end
        end
    end
    finish!(p)

    # Merge thread-local results
    for i in 1:ngp
        min_dist, min_idx = findmin([abs(dist_local[tid][i]) for tid in 1:nthreads])
        dist[i] = min_dist
        xp[:, i] = xp_local[min_idx][:, i]
    end
    
    # Apply distance truncation if needed
    for i in eachindex(dist)
        if abs(dist[i]) > norm(grid.cell_size)
            dist[i] = sign(dist[i]) * norm(grid.cell_size)
        end
    end
    
    # INFO: Analyze and export results:
    #
    # Xg, Xp, mean_PD, max_PD = SelectProjectedNodes(mesh, grid, xp, points)
    # println("Mean of projected distance: ", mean_PD)
    # println("Maximum projected distance: ", max_PD)
    #
    # # Export visualization files
    # if !isempty(Xg) && !isempty(Xp)
    #     nnp = size(Xg, 1)
    #
    #     IEN = [[i; i + nnp] for i = 1:nnp]
    #     X = vcat(Xg, Xp)
    #
    #     stl2sdf.exportToVTU("lines_STL.vtu", X, IEN, 3)
    #
    #     IEN = [[i] for i = 1:nnp]
    #     stl2sdf.exportToVTU("Xg_STL.vtu", Xg, IEN, 1)
    #     stl2sdf.exportToVTU("Xp_STL.vtu", Xp, IEN, 1)
    # end
    
    return dist, xp
end

# Helper struct to store preprocessed triangle data
struct TriangleData
    vertices::Vector{Vector{Float64}}
    edges::Vector{Vector{Float64}}
    edge_lengths::Vector{Float64}
    normal::Vector{Float64}
    center::Vector{Float64}
    aabb_min::Vector{Float64}
    aabb_max::Vector{Float64}
end

# Preprocess triangle mesh to extract and cache useful properties
function prepare_triangle_data(mesh::TriangularMesh)
    triangle_data = Vector{TriangleData}(undef, mesh.nel)
    
    for el in 1:mesh.nel
        # Extract triangle vertices
        vtx_indices = mesh.IEN[:, el]
        vertices = [mesh.X[:, vtx_indices[1]], mesh.X[:, vtx_indices[2]], mesh.X[:, vtx_indices[3]]]
        
        # Compute triangle edges
        edges = [
            vertices[2] - vertices[1],
            vertices[3] - vertices[2],
            vertices[1] - vertices[3]
        ]
        
        # Compute edge lengths
        edge_lengths = [norm(edge) for edge in edges]
        
        # Compute normal
        normal = cross(edges[1], -edges[3])
        normal = normal / norm(normal)
        
        # Compute center point
        center = (vertices[1] + vertices[2] + vertices[3]) / 3.0
        
        # Compute axis-aligned bounding box
        aabb_min = [minimum([v[i] for v in vertices]) for i in 1:3]
        aabb_max = [maximum([v[i] for v in vertices]) for i in 1:3]
        
        # Store processed data
        triangle_data[el] = TriangleData(
            vertices,
            edges,
            edge_lengths,
            normal,
            center,
            aabb_min,
            aabb_max
        )
    end
    
    return triangle_data
end

# Efficient check if a point is within range of an AABB
function point_in_box_range(point, box_min, box_max, max_dist)
    for i in 1:length(point)
        if point[i] < box_min[i] - max_dist || point[i] > box_max[i] + max_dist
            return false
        end
    end
    return true
end

 
