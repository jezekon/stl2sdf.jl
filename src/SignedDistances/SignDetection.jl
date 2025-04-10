"""
    SignDetection(mesh::Mesh, grid::Grid, points::Matrix) -> Vector{Float64}

Determine whether each point in a grid is inside or outside a tetrahedral mesh.

This function assigns a sign (+1 or -1) to each grid point based on whether it is 
inside (+1) or outside (-1) the tetrahedral mesh. It uses spatial acceleration 
structures for efficiency and parallel processing for speed.

# Arguments
- `mesh::Mesh`: The tetrahedral mesh to check against
- `grid::Grid`: The grid structure defining the domain
- `points::Matrix`: Matrix of point coordinates (3×n format)

# Returns
- `Vector{Float64}`: Vector of signs (+1 for inside, -1 for outside)
"""
function SignDetection(mesh::Mesh, grid::Grid, points::Matrix)
    X = mesh.X
    IEN = mesh.IEN
    nel = mesh.nel
    ngp = grid.ngp
    signs = -1.0 * ones(ngp)
    
    # Extract grid dimensions and properties
    grid_dims = grid.N .+ 1  # Number of points in each dimension
    grid_min = grid.AABB_min
    cell_size = grid.cell_size
    
    # Create a grid-to-tetrahedra mapping
    println("Building spatial acceleration structure...")
    grid_tetrahedra = create_grid_tetrahedra_mapping(mesh, grid, grid_dims)
    
    # Process grid points in parallel with better locality
    p_nodes = Progress(ngp, 1, "Computing signs: ", 30)
    counter_nodes = Atomic{Int}(0)
    update_interval_nodes = max(1, div(ngp, 100))
    
    # Process points in small batches to improve cache coherence
    batch_size = 64
    num_batches = ceil(Int, ngp / batch_size)
    
    @threads for batch in 1:num_batches
        start_idx = (batch - 1) * batch_size + 1
        end_idx = min(batch * batch_size, ngp)
        
        for i in start_idx:end_idx
            # Convert linear index to 3D grid index
            x = @view points[:, i]
            grid_idx = point_to_grid_index(x, grid_min, cell_size, grid_dims)
            
            # Skip invalid indices (points outside the grid bounds)
            if any(grid_idx .< 1) || any(grid_idx .> grid_dims)
                continue
            end
            
            # Get candidate tetrahedra for this grid cell
            cell_tetrahedra = grid_tetrahedra[grid_idx...]
            
            # Check if point is inside any candidate tetrahedron
            for el in cell_tetrahedra
                tetrahedron = @view X[:, IEN[:, el]]
                
                if is_point_in_tetrahedron_fastest(tetrahedron, x)
                    signs[i] = 1.0
                    break
                end
            end
            
            # Update progress less frequently to reduce atomic operation overhead
            if i == end_idx
                count = atomic_add!(counter_nodes, end_idx - start_idx + 1)
                if Threads.threadid() == 1 && (count % update_interval_nodes <= batch_size)
                    update!(p_nodes, min(count, ngp))
                end
            end
        end
    end
    
    finish!(p_nodes)
    return signs
end

"""
    create_grid_tetrahedra_mapping(mesh::Mesh, grid::Grid, grid_dims) -> Array{Vector{Int}}

Create a spatial acceleration structure that maps grid cells to potentially overlapping tetrahedra.

This function builds an efficient spatial mapping that allows quick identification of which
tetrahedra might contain a given point, greatly reducing the number of point-in-tetrahedron
tests needed.

# Arguments
- `mesh::Mesh`: The tetrahedral mesh
- `grid::Grid`: The grid structure
- `grid_dims`: Dimensions of the grid

# Returns
- `Array{Vector{Int}}`: 3D array where each cell contains indices of potentially overlapping tetrahedra
"""
function create_grid_tetrahedra_mapping(mesh::Mesh, grid::Grid, grid_dims)
    X = mesh.X
    IEN = mesh.IEN
    nel = mesh.nel
    grid_min = grid.AABB_min
    cell_size = grid.cell_size
    
    # Number of threads for parallel processing
    num_threads = Threads.nthreads()
    
    # Create thread-local storage - one grid per thread
    local_grids = [
        [Vector{Int}() for _ in 1:grid_dims[1], _ in 1:grid_dims[2], _ in 1:grid_dims[3]]
        for _ in 1:num_threads
    ]
    
    @threads for el in 1:nel
        tid = Threads.threadid()
        
        # Get tetrahedron vertices
        tet_vertices = @view X[:, IEN[:, el]]
        
        # Compute AABB of the tetrahedron
        min_bounds = vec(minimum(tet_vertices, dims=2))
        max_bounds = vec(maximum(tet_vertices, dims=2))
        
        # Convert to grid cell indices (with padding for safety)
        min_idx = max.(1, floor.(Int, (min_bounds .- grid_min) ./ cell_size) .- 1)
        max_idx = min.(grid_dims, ceil.(Int, (max_bounds .- grid_min) ./ cell_size) .+ 1)
        
        # Add this tetrahedron to all overlapping grid cells in thread-local grid
        for i in min_idx[1]:max_idx[1]
            for j in min_idx[2]:max_idx[2]
                for k in min_idx[3]:max_idx[3]
                    push!(local_grids[tid][i, j, k], el)
                end
            end
        end
    end

    # Initialize final grid
    grid_tetrahedra = [Vector{Int}() for _ in 1:grid_dims[1], _ in 1:grid_dims[2], _ in 1:grid_dims[3]]
    
    # Merge results from all thread-local grids
    for i in 1:grid_dims[1]
        for j in 1:grid_dims[2]
            for k in 1:grid_dims[3]
                for tid in 1:num_threads
                    append!(grid_tetrahedra[i, j, k], local_grids[tid][i, j, k])
                end
            end
        end
    end
    
    # Calculate and print statistics for diagnostics
    cell_counts = [length(cell) for cell in grid_tetrahedra]
    max_tets = maximum(cell_counts)
    empty_cells = count(c -> c == 0, cell_counts)
    occupied_cells = length(cell_counts) - empty_cells
    
    # println("Grid-to-tetrahedra mapping statistics:")
    # println("  Total cells: $(length(grid_tetrahedra))")
    # println("  Occupied cells: $occupied_cells ($(round(100*occupied_cells/length(grid_tetrahedra), digits=2))%)")
    # println("  Max tetrahedra per cell: $max_tets")
    avg_per_cell = occupied_cells > 0 ? sum(cell_counts) / occupied_cells : 0
    # println("  Average tetrahedra per occupied cell: $(round(avg_per_cell, digits=2))")
    
    return grid_tetrahedra
end

"""
    point_to_grid_index(point, grid_min, cell_size, grid_dims) -> Vector{Int}

Convert a 3D point to grid indices.

# Arguments
- `point`: The 3D point coordinates
- `grid_min`: Minimum coordinates of the grid
- `cell_size`: Size of grid cells (scalar or vector)
- `grid_dims`: Dimensions of the grid

# Returns
- `Vector{Int}`: 3D index of the grid cell containing the point
"""
function point_to_grid_index(point, grid_min, cell_size, grid_dims)
    # Calculate cell indices
    if isa(cell_size, Number)
        idx = floor.(Int, (point .- grid_min) ./ cell_size) .+ 1
    else
        idx = floor.(Int, (point .- grid_min) ./ cell_size) .+ 1
    end
    
    # Ensure indices are within bounds
    return max.(1, min.(grid_dims, idx))
end

"""
    is_point_in_tetrahedron_fastest(tetrahedron, point, tolerance=1e-10) -> Bool

Efficiently check if a point is inside a tetrahedron.

This optimized function uses a combination of AABB testing and barycentric coordinates
to determine if a point lies inside a tetrahedron.

# Arguments
- `tetrahedron`: Matrix of tetrahedron vertex coordinates (3×4)
- `point`: The point to test (3D vector)
- `tolerance`: Numerical tolerance for boundary cases

# Returns
- `Bool`: True if the point is inside the tetrahedron, false otherwise
"""
function is_point_in_tetrahedron_fastest(tetrahedron::AbstractMatrix{Float64}, point::AbstractVector{Float64}, tolerance::Float64=1e-10)
    # Quick AABB test
    min_bounds = minimum(tetrahedron, dims=2)
    max_bounds = maximum(tetrahedron, dims=2)
    
    if any(point .< min_bounds .- tolerance) || any(point .> max_bounds .+ tolerance)
        return false
    end
    
    # Get vertices
    v1 = tetrahedron[:, 1]
    v2 = tetrahedron[:, 2]
    v3 = tetrahedron[:, 3] 
    v4 = tetrahedron[:, 4]
    
    # Compute barycentric coordinates
    T = [v1 v2 v3 v4; 1 1 1 1]
    b = [point; 1.0]
    
    # Solve the linear system
    λ = T \ b
    
    # Point is inside if all barycentric coordinates are in [0,1] with tolerance
    return all(λ .>= -tolerance) && all(λ .<= 1.0 + tolerance)
end
