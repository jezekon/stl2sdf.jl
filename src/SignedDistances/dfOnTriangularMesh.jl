"""
    evalDistancesOnTriMesh(mesh::TriangularMesh, grid::Grid, points::Matrix{Float64})

Compute the unsigned distance field from a regular grid to a triangular mesh surface.

This function efficiently computes the minimum distance from each grid point to the
triangular mesh surface using spatial acceleration structures (KD-trees) to avoid
checking all triangles for each point.

# Arguments
- `mesh::TriangularMesh`: The triangular mesh representing the surface
- `grid::Grid`: The regular grid structure defining the domain
- `points::Matrix{Float64}`: Matrix of grid point coordinates (3×ngp)

# Returns
- `dist::Vector{Float64}`: Vector of unsigned distances from each grid point to the mesh surface
- `xp::Matrix{Float64}`: Matrix of closest point projections on the mesh for each grid point (3×ngp)
"""
function evalDistancesOnTriMesh(mesh::TriangularMesh, grid::Grid, points::Matrix{Float64})
    println("Building accelerated spatial structure...")
    # Convert the triangular mesh into a format suitable for KD-tree operations
    # Preprocessing all triangles to create an efficient searchable structure
    triangle_data = prepare_triangle_data(mesh)

    # println("Creating KD-tree for mesh vertices...")
    # Create a KD-tree from all mesh vertices for nearest neighbor searches
    vtx_kdtree = KDTree(mesh.X)

    # Convert the triangle data to a format suitable for range searches
    triangle_centers =
        Matrix{Float64}(hcat([triangle_data[i].center for i = 1:length(triangle_data)]...))

    # println("Creating KD-tree for triangle centers...")
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
    dist_local = [fill(big, ngp) for _ = 1:nthreads]
    xp_local = [zeros(Float64, nsd, ngp) for _ = 1:nthreads]

    # Create a progress bar
    p = Progress(ngp, 1, "Computing distances: ", 30)

    # Use atomic counter for progress tracking
    counter = Atomic{Int}(0)
    update_interval = max(1, div(ngp, 100))

    # Process points in batches to improve cache efficiency
    batch_size = 128
    num_batches = ceil(Int, ngp / batch_size)

    # Parallel processing of points
    Threads.@threads for batch = 1:num_batches
        tid = Threads.threadid()
        start_idx = (batch - 1) * batch_size + 1
        end_idx = min(batch * batch_size, ngp)

        for v = start_idx:end_idx
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
                if !point_in_box_range(
                    point,
                    tri_data.aabb_min,
                    tri_data.aabb_max,
                    min_dist,
                )
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
    for i = 1:ngp
        min_dist, min_idx = findmin([abs(dist_local[tid][i]) for tid = 1:nthreads])
        dist[i] = min_dist
        xp[:, i] = xp_local[min_idx][:, i]
    end

    # Apply distance truncation if needed
    for i in eachindex(dist)
        if abs(dist[i]) > norm(grid.cell_size)
            dist[i] = sign(dist[i]) * norm(grid.cell_size)
        end
    end

    return dist, xp
end

"""
    barycentricCoordinates(x₁, x₂, x₃, n, x)

Compute barycentric coordinates of point x with respect to triangle with vertices x₁, x₂, x₃.

# Arguments
- `x₁::Vector{Float64}`: First vertex of the triangle
- `x₂::Vector{Float64}`: Second vertex of the triangle
- `x₃::Vector{Float64}`: Third vertex of the triangle
- `n::Vector{Float64}`: Unit normal to the triangle face
- `x::Vector{Float64}`: Point for which to compute barycentric coordinates

# Returns
- `λ::Vector{Float64}`: Barycentric coordinates [λ₁, λ₂, λ₃]
"""
function barycentricCoordinates(
    x₁::Vector{Float64}, # First vertex of the triangle
    x₂::Vector{Float64}, # Second vertex of the triangle
    x₃::Vector{Float64}, # Third vertex of the triangle
    n::Vector{Float64},  # Unit normal to the triangle face
    x,
)                  # Point for which to compute barycentric coordinates

    # Matrix for solving the barycentric coordinates
    A = [
        (x₁[2]*n[3]-x₁[3]*n[2]) (x₂[2]*n[3]-x₂[3]*n[2]) (x₃[2]*n[3]-x₃[3]*n[2])
        (x₁[3]*n[1]-x₁[1]*n[3]) (x₂[3]*n[1]-x₂[1]*n[3]) (x₃[3]*n[1]-x₃[1]*n[3])
        (x₁[1]*n[2]-x₁[2]*n[1]) (x₂[1]*n[2]-x₂[2]*n[1]) (x₃[1]*n[2]-x₃[2]*n[1])
    ]

    # Right-hand side vector
    b = [x[2] * n[3] - x[3] * n[2], x[3] * n[1] - x[1] * n[3], x[1] * n[2] - x[2] * n[1]]

    # Find the component with maximum normal value for numerical stability
    n_max, i_max = findmax(abs.(n))
    A[i_max, :] = [1.0 1.0 1.0]
    b[i_max] = 1.0

    # Solve the system to get barycentric coordinates
    return λ = A \ b
end

"""
    SelectProjectedNodes(mesh, grid, xp, points)

Select the grid points that have been successfully projected onto the mesh surface.
This is useful for post-processing and visualization of the projection results.

# Arguments
- `mesh::TriangularMesh`: The triangular mesh
- `grid::Grid`: The regular grid structure
- `xp::Matrix{Float64}`: Matrix of projected points on the mesh (3×ngp)
- `points::Matrix{Float64}`: Matrix of original grid points (3×ngp)

# Returns
- `X::Vector{Vector{Float64}}`: Vector of grid points that have been projected
- `Xp::Vector{Vector{Float64}}`: Vector of corresponding projected points on the mesh
- `mean_PD::Float64`: Mean projection distance
- `max_PD::Float64`: Maximum projection distance
"""
function SelectProjectedNodes(
    mesh::TriangularMesh,
    grid::Grid,
    xp::Matrix{Float64},
    points::Matrix{Float64},
)

    ngp = grid.ngp # number of nodes in grid
    nsd = mesh.nsd # number of spatial dimensions

    # Initialize arrays to store valid projection pairs
    X = Vector{Vector{Float64}}()
    Xp = Vector{Vector{Float64}}()

    # Select points with valid projections (non-zero)
    for i = 1:ngp
        if sum(abs.(xp[:, i])) > 1.0e-10
            push!(X, points[:, i])
            push!(Xp, xp[:, i])
        end
    end

    # If no points were projected, return empty results with NaN statistics
    if isempty(X)
        return X, Xp, NaN, NaN
    end

    # Calculate projection distance statistics
    distances = [norm(X[i] - Xp[i]) for i = 1:length(X)]
    mean_PD = mean(distances)
    max_PD = maximum(distances)

    return X, Xp, mean_PD, max_PD
end

# Helper struct to store preprocessed triangle data for efficient distance calculations
struct TriangleData
    vertices::Vector{Vector{Float64}}  # Triangle vertices [v₁, v₂, v₃]
    edges::Vector{Vector{Float64}}     # Edge vectors [e₁, e₂, e₃]
    edge_lengths::Vector{Float64}      # Length of each edge
    normal::Vector{Float64}            # Unit normal to the triangle face
    center::Vector{Float64}            # Center point of the triangle
    aabb_min::Vector{Float64}          # Minimum corner of axis-aligned bounding box
    aabb_max::Vector{Float64}          # Maximum corner of axis-aligned bounding box
end

"""
    prepare_triangle_data(mesh::TriangularMesh)

Preprocess the triangular mesh to extract and cache useful properties for each triangle.
This significantly speeds up distance calculations by avoiding redundant computations.

# Arguments
- `mesh::TriangularMesh`: The triangular mesh

# Returns
- `Vector{TriangleData}`: Vector of preprocessed triangle data
"""
function prepare_triangle_data(mesh::TriangularMesh)
    triangle_data = Vector{TriangleData}(undef, mesh.nel)

    for el = 1:mesh.nel
        # Extract triangle vertices
        vtx_indices = mesh.IEN[:, el]
        vertices = [
            mesh.X[:, vtx_indices[1]],
            mesh.X[:, vtx_indices[2]],
            mesh.X[:, vtx_indices[3]],
        ]

        # Compute triangle edges
        edges = [
            vertices[2] - vertices[1],
            vertices[3] - vertices[2],
            vertices[1] - vertices[3],
        ]

        # Compute edge lengths
        edge_lengths = [norm(edge) for edge in edges]

        # Compute normal
        normal = cross(edges[1], -edges[3])
        normal = normal / norm(normal)

        # Compute center point
        center = (vertices[1] + vertices[2] + vertices[3]) / 3.0

        # Compute axis-aligned bounding box
        aabb_min = [minimum([v[i] for v in vertices]) for i = 1:3]
        aabb_max = [maximum([v[i] for v in vertices]) for i = 1:3]

        # Store processed data
        triangle_data[el] =
            TriangleData(vertices, edges, edge_lengths, normal, center, aabb_min, aabb_max)
    end

    return triangle_data
end

"""
    point_in_box_range(point, box_min, box_max, max_dist)

Check if a point is within a specified distance of an axis-aligned bounding box (AABB).
This is used for quick rejection testing before more expensive distance calculations.

# Arguments
- `point::Vector{Float64}`: The point to test
- `box_min::Vector{Float64}`: Minimum corner of the AABB
- `box_max::Vector{Float64}`: Maximum corner of the AABB
- `max_dist::Float64`: Maximum distance to consider

# Returns
- `Bool`: True if the point is within range of the AABB, false otherwise
"""
function point_in_box_range(point, box_min, box_max, max_dist)
    for i = 1:length(point)
        if point[i] < box_min[i] - max_dist || point[i] > box_max[i] + max_dist
            return false
        end
    end
    return true
end
