function compute_barycentric_coordinates(tetrahedron::Matrix{Float64}, point::AbstractVector{Float64})
    # Create matrix of tetrahedron vertices and point
    # For a tetrahedron with vertices v₁, v₂, v₃, v₄
    # and a point p, we solve the system:
    # p = λ₁v₁ + λ₂v₂ + λ₃v₃ + λ₄v₄
    # where λ₁ + λ₂ + λ₃ + λ₄ = 1
    
    # Matrix form: [v₁ v₂ v₃ v₄][λ₁ λ₂ λ₃ λ₄]ᵀ = [p 1]ᵀ
    # Augment vertices with a row of ones
    A = [tetrahedron; ones(1, 4)]
    b = [point; 1.0]
    
    # Solve for barycentric coordinates
    return A \ b
end

function is_point_in_tetrahedron(tetrahedron::Matrix{Float64}, point::AbstractVector{Float64}, tolerance::Float64=1e-10)
    # Compute barycentric coordinates
    bary_coords = compute_barycentric_coordinates(tetrahedron, point)
    
    # Check if all coordinates are in [0,1] with some tolerance
    return all(bary_coords .>= -tolerance) && all(bary_coords .<= 1.0 + tolerance)
end

# Function to create AABB from a set of points
function compute_aabb(points::SubArray)
  min_bounds = minimum(points, dims=2)
  max_bounds = maximum(points, dims=2)
  return min_bounds, max_bounds
end

# Function to check if a point is inside the AABB
function is_point_inside_aabb(x::SubArray, min_bounds, max_bounds)
  return all(min_bounds .<= x) && all(x .<= max_bounds)
end

# Main function for sign detection:
function SignDetection(mesh::Mesh, grid::Grid, points::Matrix)
    X = mesh.X
    IEN = mesh.IEN
    nel = mesh.nel
    ngp = grid.ngp
    signs = -1.0 * ones(ngp)
    
    # Pre-compute AABBs for all elements
    element_aabbs = Vector{NTuple{2,Vector{Float64}}}(undef, nel)
    @threads for el in 1:nel
        aabb = compute_aabb(@view X[:, IEN[:, el]])
        element_aabbs[el] = (vec(aabb[1]), vec(aabb[2]))
    end
    
    p_nodes = Progress(ngp, 1, "Processing grid nodes: ", 30)
    counter_nodes = Atomic{Int}(0)
    update_interval_nodes = max(1, div(ngp, 100))
    
    @threads for i in 1:ngp
        x = @view points[:, i]
        
        # Find potential elements that contain the point using AABB check
        candidate_elements = [el for el in 1:nel if is_point_inside_aabb(x, element_aabbs[el]...)]
        
        # Skip if no candidate elements
        if isempty(candidate_elements)
            continue
        end
        
        # Check each candidate element
        for el in candidate_elements
            tetrahedron = X[:, IEN[:, el]]
            
            # Check if point is inside tetrahedron using barycentric coordinates
            if is_point_in_tetrahedron(tetrahedron, x)
                signs[i] = 1.0
                break  # Exit element loop once we find a containing element
            end
        end
        
        # Update progress bar
        count = atomic_add!(counter_nodes, 1)
        if count % update_interval_nodes == 0 && Threads.threadid() == 1
            update!(p_nodes, count)
        end
    end
    
    finish!(p_nodes)
    return signs
end
