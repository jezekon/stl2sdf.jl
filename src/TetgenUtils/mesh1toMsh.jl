"""
    tetgen_to_msh(base_filename::String, output_filename::String)

Convert TetGen mesh files (.node and .ele) to Gmsh format (.msh).
- `base_filename`: Base name of TetGen files without extension
- `output_filename`: Output filename without extension
"""
function tetgen_to_msh(base_filename::String, output_filename::String)
    # ===== Reading nodes (.node file) =====
    node_file = open("$base_filename.node", "r")
    lines = readlines(node_file)
    close(node_file)

    # Process header to get number of points
    header = split(lines[1])
    num_points::Int = parse(Int, header[1])

    # Allocate array for points
    points = zeros(Float64, 3, num_points)

    # Determine if indexing starts from 0 or 1 based on first point
    first_point_line = split(lines[2])
    first_index = parse(Int, first_point_line[1])
    index_offset = (first_index == 0) ? 1 : 0  # If indices start from 0, add offset of 1

    # Load points
    for i = 1:num_points
        values = split(lines[i+1])
        # Convert TetGen index to array index (0-based -> 1-based if needed)
        point_idx = parse(Int, values[1]) + index_offset
        if point_idx <= num_points
            points[1, point_idx] = parse(Float64, values[2])  # x-coordinate
            points[2, point_idx] = parse(Float64, values[3])  # y-coordinate
            points[3, point_idx] = parse(Float64, values[4])  # z-coordinate
        else
            println("Warning: Point index $point_idx is out of range (max $num_points)")
        end
    end

    # ===== Reading elements (.ele file) =====
    ele_file = open("$base_filename.ele", "r")
    lines = readlines(ele_file)
    close(ele_file)

    # Process header to get number of tetrahedra
    header = split(lines[1])
    num_tets::Int = parse(Int, header[1])
    nodes_per_tet = length(split(lines[2])) - 1  # Number of nodes per tetrahedron

    # Load tetrahedra
    tetrahedra = Vector{Vector{Int}}(undef, num_tets)
    for i = 1:num_tets
        values = split(lines[i+1])
        # Skip element index, load vertex indices of tetrahedron
        tet_indices = [parse(Int, values[j]) + index_offset for j = 2:(nodes_per_tet+1)]
        tetrahedra[i] = tet_indices
    end

    # ===== Write Gmsh .msh file =====
    open("$output_filename.msh", "w") do file
        # Write MSH file header (version 4.1)
        write(file, "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n")

        # Write nodes section
        write(file, "\$Nodes\n")
        write(file, "1 $num_points $num_points 0\n")  # 1 entity, num_points nodes
        write(file, "3 1 0 $num_points\n")  # 3D entity, tag 1, num_points nodes

        # Write node indices
        for i = 1:num_points
            write(file, "$i\n")
        end

        # Write node coordinates
        for i = 1:num_points
            write(file, "$(points[1,i]) $(points[2,i]) $(points[3,i])\n")
        end
        write(file, "\$EndNodes\n")

        # Write elements section
        write(file, "\$Elements\n")
        write(file, "1 $num_tets $num_tets 0\n")  # 1 entity, num_tets elements
        write(file, "3 1 4 $num_tets\n")  # 3D entity, tag 1, type 4 (tetrahedron), num_tets elements

        # Write tetrahedra
        for i = 1:num_tets
            tet = tetrahedra[i]
            # Element number followed by its 4 nodes
            write(file, "$i $(tet[1]) $(tet[2]) $(tet[3]) $(tet[4])\n")
        end
        write(file, "\$EndElements\n")
    end

    println("Conversion completed. File saved as $(output_filename).msh")
end

# Usage for your specific files
base_filename = "cantilever_beam_density_iso-tri_mesh.1"
output_filename = "cantilever_beam_volume_mesh"
tetgen_to_msh(base_filename, output_filename)
