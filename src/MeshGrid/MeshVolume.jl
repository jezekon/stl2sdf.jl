"""
    calculate_tetrahedral_mesh_volume(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})

Calculate volume of a mesh composed of tetrahedral elements.

Inputs:
- X: Vector of node coordinate vectors, where each inner vector contains [x,y,z] coordinates
- IEN: Vector of element connectivity vectors, where each inner vector contains the indices of the 4 nodes
  of a tetrahedral element

Returns:
- total volume of the mesh (Float64)
"""
function calculate_mesh_volume(
    X::Vector{Vector{Float64}},
    IEN::Vector{Vector{Int64}}
)
    print_info("Computing tetrahedral mesh volume...")
    
    # Atomic variable for thread-safe volume accumulation
    total_volume = Atomic{Float64}(0.0)
    
    # Parallel processing of elements
    @threads for elem in 1:length(IEN)
        # For a tetrahedron, we need 4 nodes
        if length(IEN[elem]) != 4
            print_warning("Element $elem does not have 4 nodes. Skipping...")
            continue
        end
        
        # Extract coordinates of the 4 vertices of the tetrahedron
        v1 = X[IEN[elem][1]]
        v2 = X[IEN[elem][2]]
        v3 = X[IEN[elem][3]]
        v4 = X[IEN[elem][4]]
        
        # Calculate element volume using the scalar triple product formula
        # V = (1/6) * |((v2-v1) × (v3-v1))·(v4-v1)|
        edge1 = v2 - v1
        edge2 = v3 - v1
        edge3 = v4 - v1
        
        # Calculate volume using the scalar triple product
        # Convert vectors to SVector for better performance
        e1 = SVector{3, Float64}(edge1)
        e2 = SVector{3, Float64}(edge2)
        e3 = SVector{3, Float64}(edge3)
        
        # Volume = (1/6) * |determinant of the matrix formed by the 3 edges|
        elem_volume = abs(dot(cross(e1, e2), e3)) / 6.0
        
        # Thread-safe update of total volume
        atomic_add!(total_volume, elem_volume)
    end
    
    # Get final total volume
    final_volume = total_volume[]
    
    # Print results
    println("Total mesh volume: ", round(final_volume, sigdigits=6))
    
    return final_volume
end
