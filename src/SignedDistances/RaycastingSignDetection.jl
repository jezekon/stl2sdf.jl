"""
    RaycastingSignDetection.jl
    
Optimized implementation using ImplicitBVH.jl for acceleration.
Uses adaptive ray casting with early termination.
"""

"""
    RaycastResult

Structure to hold results from ray casting in a single direction.
"""
struct RaycastResult
    intersections::Int
    direction::Vector{Float64}
    is_inside::Bool
end

"""
    TriangleBVH

Wrapper structure holding BVH and original triangles for ray casting.
"""
struct TriangleBVH
    bvh::BVH
    triangles::Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
    bounding_spheres::Vector{BSphere{Float64}}
end

"""
    build_triangle_bvh(triangles)

Build BVH structure from triangles for accelerated ray casting.
"""
function build_triangle_bvh(
    triangles::Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}},
)::TriangleBVH
    n_triangles = length(triangles)
    bounding_spheres = Vector{BSphere{Float64}}(undef, n_triangles)

    # Create bounding sphere for each triangle
    @threads for i = 1:n_triangles
        v1, v2, v3 = triangles[i]
        center = (v1 + v2 + v3) / 3.0

        # Calculate radius to enclose triangle
        radius = max(norm(v1 - center), norm(v2 - center), norm(v3 - center))
        # Add small padding for numerical robustness
        bounding_spheres[i] = BSphere(center, radius * 1.01)
    end

    # Build BVH with BBox for internal nodes (more efficient for rays)
    bvh = BVH(bounding_spheres, BBox{Float64}, UInt32)

    return TriangleBVH(bvh, triangles, bounding_spheres)
end

"""
    ray_triangle_intersection(origin, direction, v1, v2, v3, ε=1e-10)

Compute ray-triangle intersection using Möller-Trumbore algorithm.
Returns parameter t if intersection exists, nothing otherwise.
"""
function ray_triangle_intersection(
    origin::Vector{Float64},
    direction::Vector{Float64},
    v1::Vector{Float64},
    v2::Vector{Float64},
    v3::Vector{Float64},
    ε::Float64 = 1e-10,
)::Union{Float64,Nothing}
    # Edge vectors
    edge1 = v2 - v1
    edge2 = v3 - v1

    # Begin calculating determinant
    h = cross(direction, edge2)
    a = dot(edge1, h)

    # Ray is parallel to triangle
    if abs(a) < ε
        return nothing
    end

    f = 1.0 / a
    s = origin - v1
    u = f * dot(s, h)

    # Check if intersection lies outside triangle
    if u < 0.0 || u > 1.0
        return nothing
    end

    q = cross(s, edge1)
    v = f * dot(direction, q)

    if v < 0.0 || u + v > 1.0
        return nothing
    end

    # Calculate parameter t
    t = f * dot(edge2, q)

    # Only return positive t (ray, not line)
    return t > ε ? t : nothing
end

"""
    check_ray_inside_bvh(origin, direction, tri_bvh, ε=1e-8)

Check if point is inside mesh using early termination.
Returns true after finding first valid intersection (odd parity).
"""
function check_ray_inside_bvh(
    origin::Vector{Float64},
    direction::Vector{Float64},
    tri_bvh::TriangleBVH,
    ε::Float64 = 1e-8,
)::Bool
    # Reshape for ImplicitBVH interface (3×1 matrices)
    points = reshape(origin, 3, 1)
    directions = reshape(direction, 3, 1)

    # Find potential intersections using BVH
    traversal = traverse_rays(tri_bvh.bvh, points, directions)

    # Collect valid intersections
    intersections = Float64[]

    for (tri_idx, _) in traversal.contacts
        v1, v2, v3 = tri_bvh.triangles[tri_idx]
        t = ray_triangle_intersection(origin, direction, v1, v2, v3, ε)
        if t !== nothing && t > ε
            push!(intersections, t)
        end
    end

    if isempty(intersections)
        return false
    end

    # Remove duplicates for robust parity check
    sort!(intersections)
    unique_count = 1

    for i = 2:length(intersections)
        if abs(intersections[i] - intersections[i-1]) > ε
            unique_count += 1
        end
    end

    # Odd count means inside
    return unique_count % 2 == 1
end

"""
    generate_ray_directions_prioritized()

Generate ray directions in priority order for adaptive sampling.
Returns ordered list where early directions give good coverage.
"""
function generate_ray_directions_prioritized()::Vector{Vector{Float64}}
    directions = Vector{Vector{Float64}}()

    # Level 1: Primary axes (6 directions)
    primary_axes = [
        [1.0, 0.0, 0.0],  # +X
        [0.0, 1.0, 0.0],  # +Y  
        [0.0, 0.0, 1.0],  # +Z
        [-1.0, 0.0, 0.0], # -X
        [0.0, -1.0, 0.0], # -Y
        [0.0, 0.0, -1.0],  # -Z
    ]

    # Level 2: Face diagonals (12 directions)
    face_diagonals = [
        [1.0, 1.0, 0.0],
        [1.0, -1.0, 0.0],
        [-1.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 0.0, -1.0],
        [-1.0, 0.0, 1.0],
        [-1.0, 0.0, -1.0],
        [0.0, 1.0, 1.0],
        [0.0, 1.0, -1.0],
        [0.0, -1.0, 1.0],
        [0.0, -1.0, -1.0],
    ]

    # Level 3: Body diagonals (8 directions)
    body_diagonals = [
        [1.0, 1.0, 1.0],
        [1.0, 1.0, -1.0],
        [1.0, -1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, -1.0, -1.0],
    ]

    # Add in priority order
    for dir in primary_axes
        push!(directions, normalize(dir))
    end

    for dir in face_diagonals
        push!(directions, normalize(dir))
    end

    for dir in body_diagonals
        push!(directions, normalize(dir))
    end

    # Additional quasi-random directions for complex cases
    # Using deterministic generation for reproducibility
    φ = (1.0 + sqrt(5.0)) / 2.0  # Golden ratio
    for i = 1:10
        θ = 2π * mod(i * φ, 1.0)
        z = 1.0 - 2.0 * (i / 10.0)
        r = sqrt(1.0 - z * z)
        push!(directions, [r * cos(θ), r * sin(θ), z])
    end

    return directions
end

"""
    adaptive_point_in_mesh_raycast(point, tri_bvh, min_rays=6, max_rays=30)

Adaptive ray casting that adjusts number of rays based on confidence.
"""
function adaptive_point_in_mesh_raycast(
    point::Vector{Float64},
    tri_bvh::TriangleBVH,
    min_rays::Int = 6,
    max_rays::Int = 30,
)::Tuple{Float64,Float64}
    directions = generate_ray_directions_prioritized()

    # Adaptive sampling
    inside_votes = 0
    total_votes = 0
    confidence_threshold = 0.85

    # Start with minimum rays
    for i = 1:min_rays
        if check_ray_inside_bvh(point, directions[i], tri_bvh)
            inside_votes += 1
        end
        total_votes += 1
    end

    # Check if we have high confidence already
    current_ratio = inside_votes / total_votes
    current_confidence = max(current_ratio, 1.0 - current_ratio)

    # Continue sampling if confidence is low
    rays_used = min_rays
    while current_confidence < confidence_threshold && rays_used < max_rays
        # Sample next batch
        batch_size = min(4, max_rays - rays_used)
        for i = (rays_used+1):(rays_used+batch_size)
            if i <= length(directions)
                if check_ray_inside_bvh(point, directions[i], tri_bvh)
                    inside_votes += 1
                end
                total_votes += 1
            end
        end
        rays_used += batch_size

        # Update confidence
        current_ratio = inside_votes / total_votes
        current_confidence = max(current_ratio, 1.0 - current_ratio)

        # Early termination for very high confidence
        if current_confidence > 0.95
            break
        end
    end

    # Determine final sign
    sign = current_ratio > 0.5 ? 1.0 : -1.0

    return (sign, current_confidence)
end

"""
    generalized_winding_number(point, triangles)

Compute generalized winding number as fallback for extremely broken meshes.
"""
function generalized_winding_number(
    point::Vector{Float64},
    triangles::Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}},
)::Float64
    winding = 0.0

    for (v1, v2, v3) in triangles
        # Vectors from point to triangle vertices
        a = normalize(v1 - point)
        b = normalize(v2 - point)
        c = normalize(v3 - point)

        # Calculate solid angle subtended by triangle
        numerator = dot(a, cross(b, c))
        denominator = 1.0 + dot(a, b) + dot(b, c) + dot(c, a)

        if abs(denominator) > 1e-10
            solid_angle = 2.0 * atan(abs(numerator), denominator)

            # Apply sign based on orientation
            if numerator < 0
                solid_angle = -solid_angle
            end

            winding += solid_angle
        end
    end

    # Normalize by 4π (total solid angle)
    winding = winding / (4π)

    # Threshold: |winding| > 0.5 means inside
    return abs(winding) > 0.5 ? 1.0 : -1.0
end

"""
    raycast_sign_detection(mesh::TriangularMesh, grid::Grid, points::Matrix{Float64}; 
                          min_rays=6, max_rays=30, confidence_threshold=0.6, 
                          use_winding_fallback=true)

Main function for optimized ray casting with adaptive sampling.
"""
function raycast_sign_detection(
    mesh::TriangularMesh,
    grid::Grid,
    points::Matrix{Float64};
    min_rays::Int = 6,
    max_rays::Int = 30,
    confidence_threshold::Float64 = 0.6,
    use_winding_fallback::Bool = true,
)::Tuple{Vector{Float64},Vector{Float64}}
    print_info("Using adaptive ray casting for sign detection...")

    # Convert mesh to triangle format
    triangles = Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}()
    for el = 1:mesh.nel
        v_indices = mesh.IEN[:, el]
        v1 = mesh.X[:, v_indices[1]]
        v2 = mesh.X[:, v_indices[2]]
        v3 = mesh.X[:, v_indices[3]]
        push!(triangles, (v1, v2, v3))
    end

    # Build BVH acceleration structure
    tri_bvh = build_triangle_bvh(triangles)

    ngp = grid.ngp
    signs = Vector{Float64}(undef, ngp)
    confidences = Vector{Float64}(undef, ngp)

    # Parallel processing
    p = Progress(ngp, 1, "Adaptive ray casting: ", 30)
    counter = Atomic{Int}(0)

    @threads for i = 1:ngp
        point = points[:, i]

        # Adaptive ray casting
        (sign, confidence) =
            adaptive_point_in_mesh_raycast(point, tri_bvh, min_rays, max_rays)

        # Fallback for low confidence
        if use_winding_fallback && confidence < confidence_threshold
            sign = generalized_winding_number(point, triangles)
            confidence = 0.8  # Higher confidence for winding number method
        end

        signs[i] = sign
        confidences[i] = confidence

        count = atomic_add!(counter, 1)
        if count % 100 == 0 && Threads.threadid() == 1
            update!(p, min(count, ngp))
        end
    end

    finish!(p)

    # Statistics
    low_confidence_count = count(c -> c < confidence_threshold, confidences)
    mean_confidence = mean(confidences)

    print_success("Ray casting completed!")
    print_info("Mean confidence: $(round(mean_confidence, digits=3))")
    print_info(
        "Low confidence points: $(low_confidence_count)/$(ngp) ($(round(100*low_confidence_count/ngp, digits=1))%)",
    )

    return (signs, confidences)
end
