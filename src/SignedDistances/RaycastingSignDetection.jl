"""
    RaycastingSignDetection.jl
    
Optimized implementation using ImplicitBVH.jl for acceleration.
Uses multi-directional ray casting with voting for robustness.
"""

using ImplicitBVH
using ImplicitBVH: BBox, BSphere
using LinearAlgebra
using Random
using Base.Threads

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
  triangles::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}
  bounding_spheres::Vector{BSphere{Float64}}
end

"""
    build_triangle_bvh(triangles)

Build BVH structure from triangles for accelerated ray casting.
"""
function build_triangle_bvh(
  triangles::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}
)::TriangleBVH
  n_triangles = length(triangles)
  bounding_spheres = Vector{BSphere{Float64}}(undef, n_triangles)

  # Create bounding sphere for each triangle
  @threads for i in 1:n_triangles
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
  ε::Float64 = 1e-10
)::Union{Float64, Nothing}
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
    count_ray_intersections_bvh(origin, direction, tri_bvh, ε=1e-8)

Count intersections using BVH acceleration. Removes duplicates.
"""
function count_ray_intersections_bvh(
  origin::Vector{Float64},
  direction::Vector{Float64},
  tri_bvh::TriangleBVH,
  ε::Float64 = 1e-8
)::Int
  # Reshape for ImplicitBVH interface (3×1 matrices)
  points = reshape(origin, 3, 1)
  directions = reshape(direction, 3, 1)

  # Find potential intersections using BVH
  traversal = traverse_rays(tri_bvh.bvh, points, directions)

  # Test actual ray-triangle intersections for candidates only
  intersections = Float64[]

  for (tri_idx, _) in traversal.contacts
    v1, v2, v3 = tri_bvh.triangles[tri_idx]
    t = ray_triangle_intersection(origin, direction, v1, v2, v3, ε)
    if t !== nothing && t > ε
      push!(intersections, t)
    end
  end

  # Remove duplicates (same distance = likely duplicate triangles)
  if length(intersections) <= 1
    return length(intersections)
  end

  sort!(intersections)
  unique_count = 1

  for i in 2:length(intersections)
    if abs(intersections[i] - intersections[i - 1]) > ε
      unique_count += 1
    end
  end

  return unique_count
end

"""
    generate_ray_directions(num_directions::Int, seed::Int=42)

Generate well-distributed ray directions for robust voting.
"""
function generate_ray_directions(
  num_directions::Int,
  seed::Int = 42
)::Vector{Vector{Float64}}
  Random.seed!(seed)

  directions = Vector{Vector{Float64}}()

  # Primary axis directions
  primary_axes = [
    [1.0, 0.0, 0.0],  # +X
    [0.0, 1.0, 0.0],  # +Y  
    [0.0, 0.0, 1.0],  # +Z
    [-1.0, 0.0, 0.0], # -X
    [0.0, -1.0, 0.0], # -Y
    [0.0, 0.0, -1.0]  # -Z
  ]

  # Diagonal directions
  diagonal_directions = [
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
    [0.0, -1.0, -1.0]
  ]

  # Add primary axes
  for dir in primary_axes
    if length(directions) < num_directions
      push!(directions, normalize(dir))
    end
  end

  # Add diagonal directions
  for dir in diagonal_directions
    if length(directions) < num_directions
      push!(directions, normalize(dir))
    end
  end

  # Fill remaining with random directions
  while length(directions) < num_directions
    θ = 2π * rand()
    φ = acos(2 * rand() - 1)

    x = sin(φ) * cos(θ)
    y = sin(φ) * sin(θ)
    z = cos(φ)

    push!(directions, normalize([x, y, z]))
  end

  return directions
end

"""
    perturb_direction(direction, perturbation_scale=1e-6)

Slightly perturb a direction to avoid degenerate cases.
"""
function perturb_direction(
  direction::Vector{Float64},
  perturbation_scale::Float64 = 1e-6
)::Vector{Float64}
  perturbation = perturbation_scale * randn(3)
  return normalize(direction + perturbation)
end

"""
    robust_point_in_mesh_raycast_bvh(point, tri_bvh, num_directions=20)

Optimized version using BVH for acceleration.
"""
function robust_point_in_mesh_raycast_bvh(
  point::Vector{Float64},
  tri_bvh::TriangleBVH,
  num_directions::Int = 20
)::Tuple{Float64, Float64}
  directions = generate_ray_directions(num_directions)

  # Process all rays in parallel
  results = Vector{RaycastResult}(undef, num_directions)

  @threads for i in 1:num_directions
    # Perturb direction slightly to avoid degeneracies
    perturbed_dir = perturb_direction(directions[i])

    # Count intersections using BVH
    intersection_count = count_ray_intersections_bvh(point, perturbed_dir, tri_bvh)

    # Odd count = inside, even count = outside
    is_inside = (intersection_count % 2 == 1)

    results[i] = RaycastResult(intersection_count, perturbed_dir, is_inside)
  end

  # Vote counting
  inside_votes = count(r -> r.is_inside, results)
  inside_fraction = inside_votes / num_directions

  # Determine sign and confidence
  sign = inside_fraction > 0.5 ? 1.0 : -1.0
  confidence = max(inside_fraction, 1.0 - inside_fraction)

  return (sign, confidence)
end

"""
    generalized_winding_number(point, triangles)

Compute generalized winding number as fallback for extremely broken meshes.
"""
function generalized_winding_number(
  point::Vector{Float64},
  triangles::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}
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
                          num_directions=20, confidence_threshold=0.6, use_winding_fallback=true)

Main function for ray casting based sign detection using BVH acceleration.
Maintains compatibility with original interface.
"""
function raycast_sign_detection(
  mesh::TriangularMesh,
  grid::Grid,
  points::Matrix{Float64};
  num_directions::Int = 20,
  confidence_threshold::Float64 = 0.6,
  use_winding_fallback::Bool = true
)::Tuple{Vector{Float64}, Vector{Float64}}
  print_info("Using BVH-accelerated ray casting for sign detection...")
  print_info("Converting mesh to triangle format...")

  # Convert mesh to triangle format
  triangles = Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
  for el in 1:mesh.nel
    v_indices = mesh.IEN[:, el]
    v1 = mesh.X[:, v_indices[1]]
    v2 = mesh.X[:, v_indices[2]]
    v3 = mesh.X[:, v_indices[3]]
    push!(triangles, (v1, v2, v3))
  end

  print_info("Building BVH structure for $(length(triangles)) triangles...")

  # Build BVH acceleration structure
  tri_bvh = build_triangle_bvh(triangles)

  ngp = grid.ngp
  signs = Vector{Float64}(undef, ngp)
  confidences = Vector{Float64}(undef, ngp)

  print_info("Processing $(ngp) grid points with $(num_directions) rays each...")

  # Progress tracking
  p = Progress(ngp, 1, "Ray casting: ", 30)
  counter = Atomic{Int}(0)
  update_interval = max(1, div(ngp, 100))

  # Process grid points in parallel
  @threads for i in 1:ngp
    point = points[:, i]

    # Primary method: BVH-accelerated ray casting with voting
    (sign, confidence) = robust_point_in_mesh_raycast_bvh(point, tri_bvh, num_directions)

    # Fallback method for low confidence points
    if use_winding_fallback && confidence < confidence_threshold
      sign = generalized_winding_number(point, triangles)
      confidence = 0.5  # Mark as fallback method
    end

    signs[i] = sign
    confidences[i] = confidence

    # Update progress
    count = atomic_add!(counter, 1)
    if count % update_interval == 0 && Threads.threadid() == 1
      update!(p, min(count, ngp))
    end
  end

  finish!(p)

  # Print statistics
  low_confidence_count = count(c -> c < confidence_threshold, confidences)
  mean_confidence = mean(confidences)

  print_success("BVH-accelerated ray casting completed!")
  print_info("Mean confidence: $(round(mean_confidence, digits=3))")
  print_info(
    "Low confidence points: $(low_confidence_count)/$(ngp) ($(round(100*low_confidence_count/ngp, digits=1))%)"
  )

  return (signs, confidences)
end
