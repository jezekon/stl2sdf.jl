"""
    run_tetgen(model_name::String, directory::String="")

Run Tetgen on an STL file to generate a tetrahedral mesh.

# Arguments
- `model_name::String`: Base name of the STL file without extension
- `directory::String="."`: Directory containing the STL file (defaults to current directory)

# Description
This function performs the following steps:
1. Clear any existing Tetgen output files for the specified model
2. Execute Tetgen with the -p option on the STL file
3. Wait for all output files to be generated (.1.node, .1.ele, .1.face, .1.edge, .1.smesh)

# Returns
- `Bool`: true if successful, false otherwise

# Example
```julia
run_tetgen("beam-approx")  # Run in current directory
run_tetgen("beam-approx", "/path/to/data")  # Run in specified directory
```
"""
function run_tetgen(model_name::String, directory::String=".")
    # Ensure directory path is correctly formatted
    dir_path = isempty(directory) ? "." : directory
    
    # Change to the directory if not the current one
    current_dir = pwd()
    try
        if dir_path != "."
            cd(dir_path)
        end
        
        # First, clean any existing Tetgen files
        clean_tetgen_files(model_name, dir_path)
        
        # Run Tetgen command
        stl_file = "$model_name.stl"
        if !isfile(stl_file)
            print_error("STL file not found: $stl_file")
            return false
        end
        
        print_info("Running Tetgen on $stl_file...")
        run_cmd = `tetgen -p $stl_file`
        
        # Execute command
        try
            run(run_cmd)
        catch e
            print_error("Error executing Tetgen: $(e)")
            return false
        end
        
        # Check for output files
        expected_extensions = [".1.node", ".1.ele", ".1.face", ".1.edge", ".1.smesh"]
        
        # Wait for files to be generated (with timeout)
        timeout = 60  # seconds
        start_time = time()
        all_files_exist = false
        
        while !all_files_exist && (time() - start_time < timeout)
            all_files_exist = true
            for ext in expected_extensions
                if !isfile("$model_name$ext")
                    all_files_exist = false
                    sleep(0.5)  # Wait before checking again
                    break
                end
            end
        end
        
        if all_files_exist
            print_success("Tetgen processing completed successfully.")
            return true
        else
            missing_files = filter(ext -> !isfile("$model_name$ext"), expected_extensions)
            print_error("Tetgen processing failed. Missing files: $(join(["$model_name$ext" for ext in missing_files], ", "))")
            return false
        end
        
    finally
        # Return to original directory
        if dir_path != "."
            cd(current_dir)
        end
    end
end

"""
    clean_tetgen_files(model_name::String, directory::String=".")

Remove any existing Tetgen output files for the specified model.

# Arguments
- `model_name::String`: Base name of the model
- `directory::String="."`: Directory containing the files
"""
function clean_tetgen_files(model_name::String, directory::String=".")
    dir_path = isempty(directory) ? "." : directory
    
    # List of possible Tetgen output extensions
    extensions = [".1.node", ".1.ele", ".1.face", ".1.edge", ".1.smesh"]
    
    for ext in extensions
        file_path = joinpath(dir_path, "$model_name$ext")
        if isfile(file_path)
            try
                rm(file_path)
                print_info("Removed existing file: $file_path")
            catch e
                print_warning("Could not remove file $file_path: $e")
            end
        end
    end
end
