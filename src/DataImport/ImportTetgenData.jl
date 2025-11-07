function import_tetgen_mesh(base_filename::String)
    # Read node file (.node)
    node_file = open("$(base_filename).node", "r")
    node_lines = readlines(node_file)
    close(node_file)

    # Parse header line
    header = split(node_lines[1])
    num_nodes = parse(Int, header[1])

    # Check if indexing starts from 0 or 1
    first_node_line = split(node_lines[2])
    first_index = parse(Int, first_node_line[1])
    index_offset = (first_index == 0) ? 1 : 0  # Add offset if zero-indexed

    # Parse node coordinates
    X = Vector{Vector{Float64}}(undef, num_nodes)
    for i = 1:num_nodes
        values = split(node_lines[i+1])
        node_idx = parse(Int, values[1]) + index_offset
        x = parse(Float64, values[2])
        y = parse(Float64, values[3])
        z = parse(Float64, values[4])
        X[node_idx] = [x, y, z]
    end

    # Read element file (.ele)
    ele_file = open("$(base_filename).ele", "r")
    ele_lines = readlines(ele_file)
    close(ele_file)

    # Parse header line
    header = split(ele_lines[1])
    num_elements = parse(Int, header[1])
    nodes_per_element = parse(Int, header[2])

    # Parse element connectivity
    IEN = Vector{Vector{Int64}}(undef, num_elements)
    for i = 1:num_elements
        values = split(ele_lines[i+1])
        element_idx = parse(Int, values[1]) + index_offset

        # Get node indices for this element
        element_nodes = Vector{Int64}(undef, nodes_per_element)
        for j = 1:nodes_per_element
            element_nodes[j] = parse(Int, values[j+1]) + index_offset
        end

        IEN[element_idx] = element_nodes
    end

    return X, IEN
end
