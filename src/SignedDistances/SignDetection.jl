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


function is_point_in_tetrahedron_halfspace(tetrahedron, point)
    # For each face of the tetrahedron
    for i in 1:4
        # Get face vertices (skipping the i-th vertex)
        vertices = [tetrahedron[:, j] for j in 1:4 if j != i]
        
        # Compute face normal pointing inward
        v1, v2, v3 = vertices
        normal = cross(v2 - v1, v3 - v1)
        
        # Ensure normal points inward (toward remaining vertex)
        remaining_vertex = tetrahedron[:, i]
        if dot(normal, remaining_vertex - v1) < 0
            normal = -normal
        end
        
        # Check if point is on the correct side of this face
        if dot(normal, point - v1) < 0
            return false  # Outside this face's half-space
        end
    end
    
    return true  # Inside all half-spaces = inside tetrahedron
end

function is_point_in_tetrahedron_volume(tetrahedron::Matrix{Float64}, point::Vector{Float64}, tolerance::Float64=1e-10)
    # Original tetrahedron vertices
    v1 = tetrahedron[:, 1]
    v2 = tetrahedron[:, 2]
    v3 = tetrahedron[:, 3]
    v4 = tetrahedron[:, 4]
    
    # Volume of original tetrahedron
    original_volume = signed_tet_volume(v1, v2, v3, v4)
    
    # For each face, create a tetrahedron with the query point
    volumes = [
        signed_tet_volume(point, v2, v3, v4),
        signed_tet_volume(v1, point, v3, v4),
        signed_tet_volume(v1, v2, point, v4),
        signed_tet_volume(v1, v2, v3, point)
    ]
    
    # If all volumes have the same sign as the original and sum to it
    # (allowing for floating point error), point is inside
    same_sign = all(sign(vol) == sign(original_volume) for vol in volumes)
    volumes_sum = sum(volumes)
    sum_approx = abs(volumes_sum - original_volume) < tolerance * abs(original_volume)
    
    return same_sign && sum_approx
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
            # if is_point_in_tetrahedron(tetrahedron, x)
            # if is_point_in_tetrahedron_halfspace(tetrahedron, x)
            if is_point_in_tetrahedron_volume(tetrahedron, x)
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
