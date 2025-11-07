"""
    import_stl(filename::String)

Import an STL file using MeshIO.jl and convert it to the expected output format.

Returns:
- X: Vector{Vector{Float64}} - A vector of vertex coordinates
- IEN: Vector{Vector{Int}} - A vector of triangle connectivity indices
"""
function import_stl(filename::String)
    try
        # Load the STL file using FileIO (which uses MeshIO under the hood)
        mesh = load(filename)

        # Extract raw vertices and faces
        raw_vertices, raw_faces = extract_raw_mesh_data(mesh)

        # Process the mesh to remove duplicates (just like the original implementation)
        @time (X, IEN) = process_stl_mesh(raw_vertices, raw_faces)

        return X, IEN
    catch e
        # Fallback to original implementation if MeshIO has issues
        @warn "MeshIO loading failed: $e"
        return fallback_import_stl(filename)
    end
end

function extract_raw_mesh_data(mesh)
    # Extract vertices 
    vertices = Vector{Vector{Float64}}()
    raw_faces = Vector{Vector{Int}}()

    # Try to extract vertices and faces based on mesh structure
    try
        # Extract vertices - converting to the expected format
        points = get_vertices(mesh)
        vertices = [Vector{Float64}([v[1], v[2], v[3]]) for v in points]

        # Extract faces - this will depend on mesh representation
        # If faces are represented as triangles with direct vertex indices
        faces = get_faces(mesh)
        for i = 1:length(faces)
            face = faces[i]
            # Create a new triangle with indices into our vertices array
            triangle = Vector{Int}([(i-1)*3 + 1, (i-1)*3 + 2, (i-1)*3 + 3])
            push!(raw_faces, triangle)
        end

        return vertices, raw_faces
    catch e
        @warn "Error extracting mesh data: $e"
        return Vector{Vector{Float64}}(), Vector{Vector{Int}}()
    end
end

function get_vertices(mesh)
    # Try different accessor methods for vertices
    for accessor in [:positions, :vertices, :points, :coordinates]
        try
            if hasproperty(mesh, accessor)
                return getproperty(mesh, accessor)
            end
        catch
            continue
        end
    end

    # Try to access through GeometryBasics methods
    try
        return coordinates(mesh)
    catch
        try
            return points(mesh)
        catch
            error("Could not extract vertices from mesh")
        end
    end
end

function get_faces(mesh)
    # Try different accessor methods for faces
    for accessor in [:faces, :elements, :triangles, :connectivity]
        try
            if hasproperty(mesh, accessor)
                return getproperty(mesh, accessor)
            end
        catch
            continue
        end
    end

    # Try to access through GeometryBasics methods
    try
        return faces(mesh)
    catch
        error("Could not extract faces from mesh")
    end
end

# Original implementation as fallback
function fallback_import_stl(filename::String)
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
        for i = 1:n_triangles
            # Skip normal vector (3 floats = 12 bytes)
            skip(file, 12)

            # Read 3 vertices (each vertex is 3 floats = 12 bytes)
            for j = 1:3
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

function process_stl_mesh(vertices, triangles)
    # Use a spatial hashing approach for much faster vertex deduplication
    tol = 1e-10
    scale_factor = floor(Int, 1/tol)  # Scaling factor for discretization

    # Dictionary to map discretized coordinates to vertex indices
    vertex_map = Dict{NTuple{3,Int},Vector{Int}}()
    unique_vertices = Vector{Vector{Float64}}()
    new_triangles = Vector{Vector{Int}}()

    # Function to discretize coordinates to integer grid for spatial hashing
    function discretize(v)
        return (
            floor(Int, v[1] * scale_factor),
            floor(Int, v[2] * scale_factor),
            floor(Int, v[3] * scale_factor),
        )
    end

    # Pre-allocate buffer for the triangle indices
    for triangle in triangles
        new_triangle = Vector{Int}(undef, 3)

        for (i, vertex_idx) in enumerate(triangle)
            vertex = vertices[vertex_idx]

            # Get discretized coordinates for spatial hashing
            grid_pos = discretize(vertex)

            # Check neighborhood cells (current cell and immediate neighbors)
            found = false
            found_idx = 0

            # Check for matching vertices in the current cell
            if haskey(vertex_map, grid_pos)
                for idx in vertex_map[grid_pos]
                    if norm(vertex - unique_vertices[idx]) < tol
                        found = true
                        found_idx = idx
                        break
                    end
                end
            end

            # If not found in current cell, check all potential neighboring cells
            if !found
                # Only check neighboring cells if vertex is close to a boundary
                for dx = -1:1, dy = -1:1, dz = -1:1
                    if dx == 0 && dy == 0 && dz == 0
                        continue  # Skip current cell, already checked
                    end

                    # Vertex position within cell
                    pos_in_cell = (
                        (vertex[1] * scale_factor) - floor(vertex[1] * scale_factor),
                        (vertex[2] * scale_factor) - floor(vertex[2] * scale_factor),
                        (vertex[3] * scale_factor) - floor(vertex[3] * scale_factor),
                    )

                    # Only check cells that could contain vertices within tolerance
                    if (
                        (dx == -1 && pos_in_cell[1] < tol * scale_factor) ||
                        (dx == 1 && pos_in_cell[1] > 1 - tol * scale_factor) ||
                        (dy == -1 && pos_in_cell[2] < tol * scale_factor) ||
                        (dy == 1 && pos_in_cell[2] > 1 - tol * scale_factor) ||
                        (dz == -1 && pos_in_cell[3] < tol * scale_factor) ||
                        (dz == 1 && pos_in_cell[3] > 1 - tol * scale_factor)
                    )

                        neighbor_pos =
                            (grid_pos[1] + dx, grid_pos[2] + dy, grid_pos[3] + dz)

                        if haskey(vertex_map, neighbor_pos)
                            for idx in vertex_map[neighbor_pos]
                                if norm(vertex - unique_vertices[idx]) < tol
                                    found = true
                                    found_idx = idx
                                    break
                                end
                            end
                        end

                        # Break out of all loops if found
                        if found
                            break
                        end
                    end
                end
            end

            # If vertex is unique, add it to our collection
            if !found
                push!(unique_vertices, vertex)
                found_idx = length(unique_vertices)

                # Add to spatial hash map
                if !haskey(vertex_map, grid_pos)
                    vertex_map[grid_pos] = [found_idx]
                else
                    push!(vertex_map[grid_pos], found_idx)
                end
            end

            # Set index in the new triangle
            new_triangle[i] = found_idx
        end

        push!(new_triangles, new_triangle)
    end

    return unique_vertices, new_triangles
end
