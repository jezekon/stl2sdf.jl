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
  # Get mesh dimensions
  nel = mesh.nel  # Number of triangles
  nen = mesh.nen  # Nodes per element (3 for triangles)
  
  # Each triangle has 3 edges
  num_edges_per_element = 3
  
  # Initialize matrix to store edge lengths
  distances = zeros(Float64, num_edges_per_element, nel)

  for e in 1:nel
    # Get vertex indices for current triangle
    v_indices = mesh.IEN[:, e]
    
    # Calculate length of each edge in the triangle
    for i in 1:num_edges_per_element
      # Get vertex indices for current edge
      start_idx = v_indices[i]
      end_idx = v_indices[mod1(i + 1, nen)]  # Wrap to first vertex after last
      
      # Get vertex coordinates
      start_point = mesh.X[:, start_idx]
      end_point = mesh.X[:, end_idx]
      
      # Store edge length
      distances[i, e] = norm(end_point - start_point)
    end
  end

  return distances
end

# Terminal color codes for highlighting important information
const RED = "\e[31m"
const BOLD = "\e[1m"
const RESET = "\e[0m"

function analyze_mesh(distances::Matrix{Float64})
  # Count number of triangular elements
  num_elements = size(distances, 2)

  # Calculate edge length statistics
  shortest_edge = minimum(distances)
  longest_edge = maximum(distances)
  avg_edge_length = mean(distances)
  median_edge_length = median(vec(distances))

  # Find smallest and largest triangles
  element_sums = sum(distances, dims=1)
  smallest_element_index = argmin(vec(element_sums))
  largest_element_index = argmax(vec(element_sums))

  # Calculate average edge length for smallest and largest triangles
  avg_edge_smallest = mean(distances[:, smallest_element_index])
  avg_edge_largest = mean(distances[:, largest_element_index])

  # Print mesh statistics with highlighting for important values
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
  # Get mesh bounding box and analyze edge lengths
  X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
  distances = calculate_edge_distances(mesh)
  analyze_mesh(distances)
  
  # Initialize grid variable
  sdf_grid = []

  # Performance note for user
  println("The time duration for 100k nodes was about 20 seconds. ")

  # Calculate number of cells based on desired cell size B
  N_new = floor(Int, maximum(X_max .- X_min) / B)
  
  # Create grid with 3-cell margin
  sdf_grid = MeshGrid.Grid(X_min, X_max, N_new, 3)
  println("Number of all grid points: ", sdf_grid.ngp)

  return sdf_grid
end

function interactive_sdf_grid_setup(mesh::TriangularMesh)
  # Get mesh bounding box and analyze edge lengths
  X_min, X_max = MeshGrid.getMesh_AABB(mesh.X)
  distances = calculate_edge_distances(mesh)
  analyze_mesh(distances)
  
  # Initialize grid variable
  sdf_grid = []

  # Performance note for user
  println("The time duration for 100k nodes was about 20 seconds. ")

  # Main interaction loop
  while true
    # Get valid grid step size from user
    while true
      print("Write a grid step based on grid analysis: ")
      user_input = readline()

      try
        B = parse(Float64, user_input)
        
        # Calculate grid dimensions based on input cell size
        N_new = floor(Int, maximum(X_max .- X_min) / B)
        sdf_grid = MeshGrid.Grid(X_min, X_max, N_new, 3)
        println("Number of all grid points: ", sdf_grid.ngp)

        break  # Exit input validation loop on success
      catch e
        if isa(e, ArgumentError)
          println("Error: Please enter a valid floating-point number.")
        else
          println("An unexpected error occurred. Please try again.")
        end
      end
    end

    # Ask if the user wants to accept current grid or try again
    print("Do you want to continue? (y/n): ")
    response = lowercase(strip(readline()))

    if response == "y"
      return sdf_grid
    else
      println("You can write a grid step.")
    end
  end
end
