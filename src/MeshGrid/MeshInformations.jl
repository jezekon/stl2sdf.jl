# Base abstract type for all mesh representations
abstract type AbstractMesh end

# Volumetric mesh representation (e.g., tetrahedral mesh)
mutable struct Mesh <: AbstractMesh
    X::Matrix{Float64}            # Matrix of node positions (dimension × number of nodes)
    IEN::Matrix{Int64}            # Element-to-node connectivity matrix
    INE::Vector{Vector{Int64}}    # Node-to-element connectivity lookup
    nsd::Int64                    # Number of spatial dimensions
    nnp::Int64                    # Total number of nodes (points)
    nen::Int64                    # Number of nodes per element
    nel::Int64                    # Total number of elements
    V_domain::Float64             # Volume of the mesh domain

    # Constructor that takes lists of nodes and element connectivity
    # and builds the complete mesh structure
    function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
        # Calculate mesh volume using tetrahedral elements
        V_domain = calculate_mesh_volume(X, IEN)

        # Convert from vector of vectors to matrices for efficient storage and access
        X = reduce(hcat, X)
        IEN = reduce(hcat, IEN)

        # Build inverse connectivity (node to elements)
        INE = nodeToElementConnectivity(X, IEN)

        # Extract mesh dimensions from the data
        nsd = size(X, 1)     # Number of dimensions from node coordinates
        nnp = size(X, 2)     # Number of nodes
        nen = size(IEN, 1)   # Nodes per element
        nel = size(IEN, 2)   # Number of elements

        return new(X, IEN, INE, nsd, nnp, nen, nel, V_domain)
    end
end

# Build node-to-element connectivity lookup for efficient neighborhood queries
function nodeToElementConnectivity(
    X::Matrix,    # Matrix of node coordinates
    IEN::Matrix,  # Element-to-node connectivity
)
    # Create array of empty vectors, one for each node
    INE = [Vector{Int64}() for _ = 1:size(X, 2)]

    # For each element and each node in that element
    for el = 1:size(IEN, 2)
        for i = 1:size(IEN, 1)
            # Add element ID to the list of elements connected to this node
            push!(INE[IEN[i, el]], el)
        end
    end
    return INE
end

# Surface triangular mesh representation for boundary representation
mutable struct TriangularMesh <: AbstractMesh
    X::Matrix{Float64}            # Matrix of node positions (dimension × number of nodes)
    IEN::Matrix{Int64}            # Element-to-node connectivity matrix
    INE::Vector{Vector{Int64}}    # Node-to-element connectivity lookup
    nsd::Int64                    # Number of spatial dimensions (3)
    nnp::Int64                    # Total number of nodes (points)
    nen::Int64                    # Number of element nodes (3)
    nel::Int64                    # Total number of elements

    # Constructor that takes lists of nodes and triangular element connectivity
    # and builds the complete triangular mesh structure
    function TriangularMesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
        # Convert from vector of vectors to matrices for efficient storage and access
        X = reduce(hcat, X)
        IEN = reduce(hcat, IEN)

        # Build inverse connectivity (node to elements)
        INE = nodeToElementConnectivity(X, IEN)

        # Extract mesh dimensions from the data
        nsd = size(X, 1)     # Number of dimensions from node coordinates
        nnp = size(X, 2)     # Number of nodes
        nen = size(IEN, 1)   # Nodes per element (3 for triangles)
        nel = size(IEN, 2)   # Number of elements

        return new(X, IEN, INE, nsd, nnp, nen, nel)
    end
end
