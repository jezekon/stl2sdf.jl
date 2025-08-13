"""
    RaycastingSignDetection.jl
    
Fallback implementation for sign detection when Tetgen fails.
Uses multi-directional ray casting with voting for robustness.
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
    ray_triangle_intersection(origin, direction, v1, v2, v3, ε=1e-10)

Compute ray-triangle intersection using Möller-Trumbore algorithm.
Returns parameter t if intersection exists, nothing otherwise.

# Arguments

  - `origin`: Ray origin point
  - `direction`: Ray direction (should be normalized)
  - `v1, v2, v3`: Triangle vertices
  - `ε`: Numerical tolerance
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
    count_ray_triangle_intersections(origin, direction, triangles, ε=1e-10)

Count intersections between a ray and collection of triangles.
Removes duplicate intersections at similar distances.

# Arguments

  - `origin`: Ray origin
  - `direction`: Ray direction (normalized)
  - `triangles`: Vector of triangles, each triangle is (v1, v2, v3)
  - `ε`: Tolerance for duplicate detection
"""
function count_ray_triangle_intersections(
  origin::Vector{Float64},
  direction::Vector{Float64},
  triangles::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}},
  ε::Float64 = 1e-8
)::Int
  intersections = Float64[]

  for (v1, v2, v3) in triangles
    t = ray_triangle_intersection(origin, direction, v1, v2, v3)
    if t !== nothing && t > ε
      push!(intersections, t)
    end
  end

  # Remove duplicates (same distance = likely duplicate triangles)
  if length(intersections) <= 1
    return length(intersections)
  end

  sort!(intersections)
  unique_intersections = [intersections[1]]

  for i in 2:length(intersections)
    if abs(intersections[i] - intersections[i - 1]) > ε
      push!(unique_intersections, intersections[i])
    end
  end

  return length(unique_intersections)
end

"""
    generate_ray_directions(num_directions::Int, seed::Int=42)

Generate well-distributed ray directions for robust voting.
Uses combination of axis-aligned, diagonal, and random directions.

# Arguments

  - `num_directions`: Total number of directions to generate
  - `seed`: Random seed for reproducibility
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
    # Generate random direction on unit sphere
    θ = 2π * rand()  # azimuthal angle
    φ = acos(2 * rand() - 1)  # polar angle

    x = sin(φ) * cos(θ)
    y = sin(φ) * sin(θ)
    z = cos(φ)

    dir = [x, y, z]
    push!(directions, normalize(dir))
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
  perturbed = direction + perturbation
  return normalize(perturbed)
end

"""
    robust_point_in_mesh_raycast(point, triangles, num_directions=20, confidence_threshold=0.6)

Determine if a point is inside a mesh using multi-directional ray casting with voting.

# Arguments

  - `point`: Point to test
  - `triangles`: Collection of triangles defining the mesh
  - `num_directions`: Number of ray directions to test
  - `confidence_threshold`: Minimum fraction of votes needed for confident decision

# Returns

  - `(sign, confidence)`: sign ∈ {-1.0, 1.0}, confidence ∈ [0, 1]
"""
function robust_point_in_mesh_raycast(
  point::Vector{Float64},
  triangles::Vector{Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}},
  num_directions::Int = 20,
  confidence_threshold::Float64 = 0.6
)::Tuple{Float64, Float64}
  directions = generate_ray_directions(num_directions)
  results = Vector{RaycastResult}()

  for dir in directions
    # Perturb direction slightly to avoid degeneracies
    perturbed_dir = perturb_direction(dir)

    # Count intersections
    intersection_count = count_ray_triangle_intersections(point, perturbed_dir, triangles)

    # Odd count = inside, even count = outside
    is_inside = (intersection_count % 2 == 1)

    push!(results, RaycastResult(intersection_count, perturbed_dir, is_inside))
  end

  # Vote counting
  inside_votes = count(r -> r.is_inside, results)
  outside_votes = length(results) - inside_votes

  total_votes = length(results)
  inside_fraction = inside_votes / total_votes

  # Determine sign
  sign = inside_fraction > 0.5 ? 1.0 : -1.0

  # Calculate confidence
  confidence = max(inside_fraction, 1.0 - inside_fraction)
  println("confidence: ", confidence)

  return (sign, confidence)
end

"""
    generalized_winding_number(point, triangles)

Compute generalized winding number as fallback for extremely broken meshes.
More robust than ray casting for non-manifold geometry.

# Arguments

  - `point`: Query point
  - `triangles`: Collection of oriented triangles

# Returns

  - `sign`: +1.0 if inside, -1.0 if outside
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

Main function for ray casting based sign detection as fallback when Tetgen fails.

# Arguments

  - `mesh`: Triangular surface mesh
  - `grid`: Grid structure
  - `points`: Matrix of grid points (3×n)
  - `num_directions`: Number of ray directions for voting
  - `confidence_threshold`: Confidence threshold for using winding number fallback
  - `use_winding_fallback`: Whether to use winding number for low-confidence points

# Returns

  - `signs`: Vector of signs (+1 inside, -1 outside)
  - `confidences`: Vector of confidence values
"""
function raycast_sign_detection(
  mesh::TriangularMesh,
  grid::Grid,
  points::Matrix{Float64};
  num_directions::Int = 20,
  confidence_threshold::Float64 = 0.6,
  use_winding_fallback::Bool = true
)::Tuple{Vector{Float64}, Vector{Float64}}
  print_info("Using ray casting fallback for sign detection...")
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

  ngp = grid.ngp
  signs = zeros(Float64, ngp)
  confidences = zeros(Float64, ngp)

  print_info("Processing $(ngp) grid points with $(num_directions) rays each...")

  # Progress tracking
  p = Progress(ngp, 1, "Ray casting: ", 30)
  counter = Atomic{Int}(0)
  update_interval = max(1, div(ngp, 100))

  # Parallel processing
  @threads for i in 1:ngp
    point = points[:, i]

    # Primary method: ray casting with voting
    (sign, confidence) =
      robust_point_in_mesh_raycast(point, triangles, num_directions, confidence_threshold)

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

  print_success("Ray casting completed!")
  print_info("Mean confidence: $(round(mean_confidence, digits=3))")
  print_info(
    "Low confidence points: $(low_confidence_count)/$(ngp) ($(round(100*low_confidence_count/ngp, digits=1))%)"
  )

  return (signs, confidences)
end
