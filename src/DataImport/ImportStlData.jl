function import_stl(filename::String)
  # Determine if the file is ASCII or binary
  is_ascii = detect_stl_type(filename)

  if is_ascii
    vertices, triangles = read_ascii_stl(filename)
  else
    vertices, triangles = read_binary_stl(filename)
  end

  # Process the vertices and triangles to create the needed data structures
  X, IEN = process_stl_mesh(vertices, triangles)

  return X, IEN
end

function detect_stl_type(filename::String)
  open(filename, "r") do file
    # Read first 5 characters
    header = read(file, 5)
    # If it starts with "solid", it's likely ASCII (though not guaranteed)
    return String(header) == "solid"
  end
end

function read_binary_stl(filename::String)
  vertices = Vector{Float64}[]
  triangles = Vector{Int64}[]

  open(filename, "r") do file
    # Skip 80-byte header
    skip(file, 80)

    # Read number of triangles (4-byte uint32)
    n_triangles = read(file, UInt32)

    # For each triangle
    for i in 1:n_triangles
      # Skip normal vector (3 floats = 12 bytes)
      skip(file, 12)

      # Read 3 vertices (each vertex is 3 floats = 12 bytes)
      for j in 1:3
        x = read(file, Float32)
        y = read(file, Float32)
        z = read(file, Float32)
        push!(vertices, [Float64(x), Float64(y), Float64(z)])
      end

      # Create triangle connectivity (using current vertices)
      base_idx = length(vertices) - 3
      push!(triangles, [base_idx + 1, base_idx + 2, base_idx + 3])

      # Skip 2-byte attribute
      skip(file, 2)
    end
  end

  return vertices, triangles
end

function read_ascii_stl(filename::String)
  vertices = Vector{Float64}[]
  triangles = Vector{Int64}[]

  open(filename, "r") do file
    lines = readlines(file)

    vertex_count = 0
    for (i, line) in enumerate(lines)
      line = strip(line)

      if startswith(line, "vertex")
        vertex_count += 1

        # Extract the vertex coordinates
        parts = split(line)
        if length(parts) >= 4
          x = parse(Float64, parts[2])
          y = parse(Float64, parts[3])
          z = parse(Float64, parts[4])
          push!(vertices, [x, y, z])

          # If we have 3 vertices, create a triangle
          if vertex_count % 3 == 0
            base_idx = length(vertices) - 2
            push!(triangles, [base_idx, base_idx + 1, base_idx + 2])
          end
        end
      end
    end
  end

  return vertices, triangles
end

# function process_stl_mesh(vertices, triangles)
#     # Process raw vertices and triangles to remove duplicates and create proper connectivity
#
#     # Create unique vertices
#     unique_vertices = Vector{Vector{Float64}}()
#     vertex_map = Dict{Vector{Float64}, Int}()
#
#     # Create new connectivity list
#     new_triangles = Vector{Vector{Int}}()
#
#     for triangle in triangles
#         new_triangle = Vector{Int}()
#
#         for vertex_idx in triangle
#             vertex = vertices[vertex_idx]
#
#             # Check if vertex already exists
#             if haskey(vertex_map, vertex)
#                 # Use existing vertex
#                 push!(new_triangle, vertex_map[vertex])
#             else
#                 # Add new unique vertex
#                 push!(unique_vertices, vertex)
#                 vertex_map[vertex] = length(unique_vertices)
#                 push!(new_triangle, length(unique_vertices))
#             end
#         end
#
#         push!(new_triangles, new_triangle)
#     end
#
#     # Format for TriangularMesh
#     X = unique_vertices
#     IEN = new_triangles
#
#     return X, IEN
# end

function process_stl_mesh(vertices, triangles)
  # Use a tolerance for vertex comparison
  tol = 1e-10

  # Create unique vertices
  unique_vertices = Vector{Vector{Float64}}()

  # Create new connectivity list
  new_triangles = Vector{Vector{Int}}()

  for triangle in triangles
    new_triangle = Vector{Int}()

    for vertex_idx in triangle
      vertex = vertices[vertex_idx]

      # Check if vertex already exists (within tolerance)
      found = false
      for (i, unique_vertex) in enumerate(unique_vertices)
        if norm(unique_vertex - vertex) < tol
          push!(new_triangle, i)
          found = true
          break
        end
      end

      if !found
        # Add new unique vertex
        push!(unique_vertices, vertex)
        push!(new_triangle, length(unique_vertices))
      end
    end

    push!(new_triangles, new_triangle)
  end

  return unique_vertices, new_triangles
end
