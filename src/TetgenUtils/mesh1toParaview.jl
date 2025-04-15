using WriteVTK
"""
    tetgen_to_vtk(base_filename::String, output_filename::String)

Convert TetGen mesh files (.node and .ele) to VTK format (.vtu).
- `base_filename`: Base name of TetGen files without extension
- `output_filename`: Output filename without extension
"""
function tetgen_to_vtk(base_filename, output_filename)
    # Loading nodes (.node file)
    node_file = open("$base_filename.node", "r")
    lines = readlines(node_file)
    close(node_file)
    
    # Processing header
    header = split(lines[1])
    num_points = parse(Int, header[1])
    
    # Allocating array for points - we need to know the exact number of points
    points = zeros(Float64, 3, num_points)
    
    # Determining if indexing starts from 0 or 1 based on the first point
    first_point_line = split(lines[2])
    first_index = parse(Int, first_point_line[1])
    index_offset = (first_index == 0) ? 1 : 0  # If indices start from 0, we add offset 1
    
    # Loading points
    for i = 1:num_points
        values = split(lines[i+1])
        # Converting index from TetGen to array index (0-based -> 1-based if needed)
        point_idx = parse(Int, values[1]) + index_offset
        if point_idx <= num_points
            points[1, point_idx] = parse(Float64, values[2])
            points[2, point_idx] = parse(Float64, values[3])
            points[3, point_idx] = parse(Float64, values[4])
        else
            println("Warning: Point index $point_idx is out of range (max $num_points)")
        end
    end
    
    # Loading elements (.ele file)
    ele_file = open("$base_filename.ele", "r")
    lines = readlines(ele_file)
    close(ele_file)
    
    # Processing header
    header = split(lines[1])
    num_tets = parse(Int, header[1])
    nodes_per_tet = length(split(lines[2])) - 1  # Number of nodes per tetrahedron
    
    # Loading tetrahedra
    cells = Array{MeshCell}(undef, num_tets)
    for i = 1:num_tets
        values = split(lines[i+1])
        # Skipping element index, loading vertex indices of tetrahedron
        tet_indices = [parse(Int, values[j]) + index_offset for j in 2:(nodes_per_tet+1)]
        cells[i] = MeshCell(VTKCellTypes.VTK_TETRA, tet_indices)
    end
    
    # Creating VTK file
    vtk_grid(output_filename, points, cells) do vtk
        # Here you can add additional attributes if available
    end
    
    println("Conversion completed. File saved as $(output_filename).vtu")
end

# Usage for your specific files
base_filename = "cantilever_beam_density_iso-tri_mesh.1"
output_filename = "cantilever_beam_volume_mesh"
tetgen_to_vtk(base_filename, output_filename)
