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
    # For each grid cell, store indices of tetrahedra that might intersect it
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
        
        local_inside_count = 0
        
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
                
                # Use the fastest point-in-tetrahedron test
                # Based on profiling, we're using is_point_in_tetrahedron_volume which
                # has shown to be the most robust despite having similar performance
                if is_point_in_tetrahedron_fastest(tetrahedron, x)
                    signs[i] = 1.0
                    local_inside_count += 1
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
    
    # Process each tetrahedron
    p_elements = Progress(nel, 1, "Mapping tetrahedra to grid: ", 30)
    
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
        
        # Update progress (only from thread 1)
        if Threads.threadid() == 1
            update!(p_elements, el)
        end
    end
    finish!(p_elements)

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
    avg_tets = mean(cell_counts)
    empty_cells = count(c -> c == 0, cell_counts)
    occupied_cells = length(cell_counts) - empty_cells
    
    println("Grid-to-tetrahedra mapping statistics:")
    println("  Total cells: $(length(grid_tetrahedra))")
    println("  Occupied cells: $occupied_cells ($(round(100*occupied_cells/length(grid_tetrahedra), digits=2))%)")
    println("  Max tetrahedra per cell: $max_tets")
    avg_per_cell = occupied_cells > 0 ? sum(cell_counts) / occupied_cells : 0
    println("  Average tetrahedra per occupied cell: $(round(avg_per_cell, digits=2))")
    
    return grid_tetrahedra
end

# Convert a 3D point to grid index
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

# Optimized point-in-tetrahedron test
function is_point_in_tetrahedron_fastest(tetrahedron::AbstractMatrix{Float64}, point::AbstractVector{Float64}, tolerance::Float64=1e-10)
    # This implementation combines efficiency and robustness
    # First, do a quick AABB test
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
    
    # Compute barycentric coordinates more efficiently
    T = [v1 v2 v3 v4; 1 1 1 1]
    b = [point; 1.0]
    
    # Use pre-computed matrix inverse if possible
    # For simplicity in this example, we'll solve directly
    λ = T \ b
    
    # Check if all barycentric coordinates are in [0,1] with tolerance
    return all(λ .>= -tolerance) && all(λ .<= 1.0 + tolerance)
end

# Efficient signed tetrahedron volume calculation
function signed_tet_volume(v1, v2, v3, v4)
    # Compute vectors from v1 to other vertices
    a = v2 - v1
    b = v3 - v1
    c = v4 - v1
    
    # Calculate volume using scalar triple product
    # Volume = (1/6) * dot(a, cross(b, c))
    return dot(a, cross(b, c)) / 6.0
end
