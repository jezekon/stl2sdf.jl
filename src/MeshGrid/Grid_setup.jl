"""
    interactive_sdf_grid_setup(mesh::Mesh)

Interactively set up a Signed Distance Function (SDF) grid based on a given mesh.

This function guides the user through the process of creating an SDF grid from a mesh,
allowing for customization of the grid step size. It performs the following steps:

1. Calculates the bounding box of the input mesh.
2. Computes and analyzes edge distances in the mesh.
3. Prompts the user to input a grid step size, with error checking for valid float input.
4. Calculates the number of grid points based on the user's input.
5. Creates an SDF grid using the specified parameters.
6. Allows the user to confirm the setup or adjust the grid step size.

The function continues to prompt for input until the user is satisfied with the grid setup.

Parameters:
- mesh::Mesh: The input mesh to be converted into an SDF grid.

Returns:
- MeshGrid.Grid: The configured SDF grid based on user input.

Note: The function informs the user that processing 100k nodes takes approximately 20 seconds,
to help in estimating computational time for larger grids.
"""

function calculate_edge_distances(mesh::TriangularMesh)
  nel = mesh.nel  # number of triangles
  nen = mesh.nen  # number of nodes per element (should be 3 for triangles)
  
  # For triangular meshes, each element has 3 edges
  num_edges_per_element = 3
  
  # Initialize output matrix: 3 edges per triangle, nel triangles
  distances = zeros(Float64, num_edges_per_element, nel)

  for e in 1:nel
    # Get vertex indices for this triangle
    v_indices = mesh.IEN[:, e]  # Should be a vector of 3 indices
    
    # For triangular mesh, iterate through its 3 edges
    for i in 1:num_edges_per_element
      # Calculate edge vertices (wrapping around for the last edge)
      start_idx = v_indices[i]
      end_idx = v_indices[mod1(i + 1, nen)]  # mod1 ensures it wraps to 1 when i = 3
      
      # Get vertex coordinates
      start_point = mesh.X[:, start_idx]
      end_point = mesh.X[:, end_idx]
      
      # Calculate distance using the norm function
      distances[i, e] = norm(end_point - start_point)
    end
  end

  return distances
end

const RED = "\e[31m"
const BOLD = "\e[1m"
const RESET = "\e[0m"

function analyze_mesh(distances::Matrix{Float64})
  num_elements = size(distances, 2)  # number of elements

  # Overall shortest and longest edge
  shortest_edge = minimum(distances)
  longest_edge = maximum(distances)

  # Average edge length
  avg_edge_length = mean(distances)

  # Median edge length
  median_edge_length = median(vec(distances))

  # Sum of edge lengths for each element
  element_sums = sum(distances, dims=1)

  # Indices of the smallest and largest elements
  smallest_element_index = argmin(vec(element_sums))
  largest_element_index = argmax(vec(element_sums))

  # Average edge length of the smallest and largest elements
  avg_edge_smallest = mean(distances[:, smallest_element_index])
  avg_edge_largest = mean(distances[:, largest_element_index])

  println("Mesh Statistics:")
  println("----------------")
  println("Number of elements: ", num_elements)
  println("Shortest edge in the mesh: $(RED)$(BOLD)$(round(shortest_edge, digits=4))$(RESET)")
  println("Longest edge in the mesh: ", round(longest_edge, digits=4))
  println("Average edge length: ", round(avg_edge_length, digits=4))
  println("Median edge length: ", round(median_edge_length, digits=4))
  println("Average edge length of the smallest element: $(RED)$(BOLD)$(round(avg_edge_smallest, digits=4))$(RESET)")
  println("Average edge length of the largest element: ", round(avg_edge_largest, digits=4))
end

function noninteractive_sdf_grid_setup(mesh::TriangularMesh, B::Float64)
  X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
  distances = calculate_edge_distances(mesh)
  analyze_mesh(distances)
  sdf_grid = []

  println("The time duration for 100k nodes was about 20 seconds. ")

  N_new = floor(Int, maximum(X_max .- X_min) / B)
  sdf_grid = MeshGrid.Grid(X_min, X_max, N_new, 3)
  println("Number of all grid points: ", sdf_grid.ngp)

  return sdf_grid
end

function interactive_sdf_grid_setup(mesh::TriangularMesh)
  X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
  distances = calculate_edge_distances(mesh)
  analyze_mesh(distances)
  sdf_grid = []

  println("The time duration for 100k nodes was about 20 seconds. ")

  while true
    while true
      print("Write a grid step based on grid analysis: ")
      user_input = readline()

      try
        B = parse(Float64, user_input)

        N_new = floor(Int, maximum(X_max .- X_min) / B)
        sdf_grid = MeshGrid.Grid(X_min, X_max, N_new, 3)
        println("Number of all grid points: ", sdf_grid.ngp)

        break  # Exit the inner loop if the input is valid
      catch e
        if isa(e, ArgumentError)
          println("Error: Please enter a valid floating-point number.")
        else
          println("An unexpected error occurred. Please try again.")
        end
      end
    end

    print("Do you want to continue? (y/n): ")
    odpoved = lowercase(strip(readline()))

    if odpoved == "y"
      return sdf_grid
    else
      println("You can write a grid step.")
    end
  end
end


