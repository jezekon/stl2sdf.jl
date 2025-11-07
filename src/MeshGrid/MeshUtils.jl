## Node position for each element
# Structure to store node coordinates for each element
# separated by x, y, z components for better cache locality
mutable struct NodalCoordinatesInElement{T<:Array}
    x::T    # x-coordinates of element nodes
    y::T    # y-coordinates of element nodes
    z::T    # z-coordinates of element nodes
end

# Extract nodal coordinates for each element in the mesh
function NodePosition3D(mesh::AbstractMesh)
    # Initialize arrays with dimensions matching element connectivity table
    table = size(mesh.IEN)
    EN_x = zeros(Float64, table) # Number of elements = number of Gauss points
    EN_y = zeros(Float64, table)
    EN_z = zeros(Float64, table)

    # For each node in each element
    for i = 1:(length(mesh.IEN[:, 1]))
        for j = 1:(length(mesh.IEN[1, :]))
            # Get node ID and its coordinates
            ID_n = Int(mesh.IEN[i, j])
            x = mesh.X[1, ID_n]
            y = mesh.X[2, ID_n]
            z = mesh.X[3, ID_n]

            # Store coordinates in respective arrays
            EN_x[i, j] = x
            EN_y[i, j] = y
            EN_z[i, j] = z
        end
    end

    # Create and return structured node position data
    EN = NodalCoordinatesInElement(EN_x, EN_y, EN_z)
    return EN # Output contains three arrays - coordinates of nodes for each element (row position = element ID)
end
