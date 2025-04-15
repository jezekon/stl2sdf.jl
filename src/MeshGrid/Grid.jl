mutable struct Grid
    AABB_min::Vector{Float64}  # Minimum coordinates of the Axis-Aligned Bounding Box (AABB)
    AABB_max::Vector{Float64}  # Maximum coordinates of the AABB
    N::Vector{Int64}           # Number of cells along the longest side (along some axis)
    cell_size::Float64         # Size of each cell in the grid
    ngp::Int64                 # Total number of grid points

    # Marginal cells are added around the AABB to provide a buffer zone.
    function Grid(
        AABB_min::Vector{Float64},
        AABB_max::Vector{Float64},
        N_max::Int64,
        margineCells::Int64 = 3,
    )
        # Calculate cell size based on the maximum dimension
        cell_size = maximum(AABB_max .- AABB_min) ./ N_max

        # Adjusting the AABB with marginal cells
        AABB_min = AABB_min .- margineCells * cell_size
        AABB_max = AABB_max .+ margineCells * cell_size

        # Calculate number of cells in each dimension
        N = ceil.(Int64, (AABB_max .- AABB_min) / cell_size)
        
        # Adjust AABB_max to ensure exact fit with cells
        AABB_max = AABB_min + N .* cell_size

        # Calculate total number of grid points
        ngp = prod(N .+ 1) # number of grid points

        return new(AABB_min, AABB_max, N, cell_size, ngp)
    end
end


# Struct for managing a linked list in the context of a grid. (Useful for spatial hashing or similar applications.)
mutable struct LinkedList # division of regular grid into regions
    grid::Grid           # The grid associated with the linked list
    head::Vector{Int64}  # Array representing the head of each list
    next::Vector{Int64}  # Array representing the next element in each list
    N::Vector{Float64}   # Number of cells along each axis of the grid

    # Constructor for LinkedList.
    # Maps points in a 3D space to the corresponding grid cells.
    function LinkedList(grid::Grid, X::Matrix{Float64})

        np = grid.ngp # Number of grid points
        N = grid.N
        AABB_min = grid.AABB_min
        AABB_max = grid.AABB_max

        # Initialize head and next arrays with -1 (indicating empty)
        head = -1 * ones(Int64, prod(N .+ 1))
        next = -1 * ones(Int64, np)

        # Calculate grid indices for each point
        I = floor.(N .* (X .- AABB_min) ./ (AABB_max .- AABB_min))
        Ia = Int.(I[3, :] .* (N[1] + 1) * (N[2] + 1) .+ I[2, :] .* (N[1] + 1) .+ I[1, :] .+ 1)

        # Construct the linked list
        for i = 1:np
            next[i] = head[Ia[i]]
            head[Ia[i]] = i
        end

        return new(grid, head, next, N)
    end
end


# Function to determine the Axis-Aligned Bounding Box (AABB) of a mesh.
function getMesh_AABB(X::Matrix{Float64})
    X_min = vec(minimum(X, dims = 2))
    X_max = vec(maximum(X, dims = 2))
    return X_min, X_max
end

# Generates grid points for a given grid structure.
# Returns a matrix where each column represents the coordinates of a grid point.
function generateGridPoints(grid::Grid)::Matrix{Float64}
    # Initialize matrix to store all grid point coordinates
    X = zeros(3, grid.ngp)
    a = 1
    # Iterate through all grid points in 3D space
    for k = 0:grid.N[3]
        for j = 0:grid.N[2]
            for i = 0:grid.N[1]
                X[:, a] = grid.AABB_min .+ grid.cell_size .* [i, j, k]
                a += 1
            end
        end
    end
    return X
end

# Generates connectivity array for hexahedral elements in the grid
function generateConnectivityArray(grid::Grid)::Vector{Vector{Int64}}
    N = grid.N .+ 1 # N is now the number of vertices in a row, not the number of cells
    # Initialize array of element connectivity
    IEN = [fill(0, 8) for _ in 1:prod(N.-1)]
    m = 1
    # Generate 8-node connectivity for each hexahedral element
    for k = 1:N[3]-1
        for j = 1:N[2]-1
            for i = 1:N[1]-1
                IEN[m] = [
                    Int((k-1) * N[1]*N[2] + (j-1) * N[1] + i ),
                    Int((k-1) * N[1]*N[2] + (j-1) * N[1] + i+1 ),
                    Int((k-1) * N[1]*N[2] + j * N[1] + i+1 ),
                    Int((k-1) * N[1]*N[2] + j * N[1] + i ),
                    Int(k * N[1]*N[2] + (j-1) * N[1] + i),
                    Int(k * N[1]*N[2] + (j-1) * N[1] + i+1 ),
                    Int(k * N[1]*N[2] + j * N[1] + i+1 ),
                    Int(k * N[1]*N[2] + j * N[1] + i ),
                ]
                m = m+1
            end
        end
    end
    #IEN = reduce(hcat, IEN) # transform Vector of Vectors to Matrix
    return IEN
end


# Calculate the mini-AABB indices within the grid for a given geometry
function calculateMiniAABB_grid(
    Xt::Matrix{Float64}, # Coordinates of the geometry (e.g., triangle)
    δ::Float64,          # Safety margin
    N::Vector{Int64},    # Grid dimensions
    AABB_min::Vector{Float64}, # Grid minimum bounds
    AABB_max::Vector{Float64}, # Grid maximum bounds
    nsd::Int)            # Number of spatial dimensions
    
    # Calculate bounds of the geometry with safety margin
    Xt_min = minimum(Xt, dims = 2) .- δ # bottom left corner coordinates
    Xt_max = maximum(Xt, dims = 2) .+ δ # top right corner coordinates

    # Mini AABB for triangle:
    # Calculate the indices in the grid for the mini-AABB
    I_min = floor.(N .* (Xt_min .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid
    I_max = floor.(N .* (Xt_max .- AABB_min) ./ (AABB_max .- AABB_min)) # Triangle location (index) within the grid

    # Ensure indices are within grid bounds
    for j = 1:nsd # check if inside grid AABB
        if (I_min[j] < 0)
            I_min[j] = 0
        end
        if (I_max[j] >= N[j])
            I_max[j] = N[j]
        end
    end

    # Create a product iterator for all grid cells that intersect with the mini-AABB
    Is = Iterators.product( # range of indices within the mini-AABB
        I_min[1]:I_max[1],
        I_min[2]:I_max[2],
        I_min[3]:I_max[3],
    )

    return Is
end
