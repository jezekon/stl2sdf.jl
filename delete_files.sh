#!/bin/bash

# This script deletes files with specific extensions (.vtu, .vtk, .jld2, .png, .vtr) in the current directory
# and in the ./test/ subdirectory if it exists. The script is designed to:
# - Search only in the current directory and ./test/ subdirectory (no deeper subdirectories)
# - Show all files that will be deleted before taking action
# - Ask for confirmation before deleting
# - Handle all three extensions (.vtu, .vtk, .vti, .jld2, .png, .vtr) in a single operation
# - Provide clear feedback about what files were found and what actions were taken
#
# Usage: ./delete_files.sh
# The script will prompt for confirmation before deleting any files.

# Function to delete files in the specified directory - using extension concatenation
delete_files() {
  local dir="$1"

  # Find all requested files at once using -o (OR) operator
  files=$(find "$dir" -maxdepth 1 -type f \( -name "*.vtu" -o -name "*.vtk" -o -name "*.vtr" -o -name "*.vti" -o -name "*.jld2" -o -name "*.png" -o -name "*.edge" -o -name "*.ele" -o -name "*.face" -o -name "*.node" -o -name "*.smesh" \))

  # If any files were found
  if [ -n "$files" ]; then
    echo "Found files to delete in directory $dir:"
    echo "$files"
    echo

    # Ask for confirmation before deletion
    read -p "Do you want to delete these files? (y/n): " answer
    if [ "$answer" = "y" ]; then
      # Delete all found files at once
      find "$dir" -maxdepth 1 -type f \( -name "*.vtu" -o -name "*.vtk" -o -name "*.vtr" -o -name "*.vti" -o -name "*.jld2" -o -name "*.png" -o -name "*.edge" -o -name "*.ele" -o -name "*.face" -o -name "*.node" -o -name "*.smesh" \) -delete
      echo "Files have been deleted."
    else
      echo "File deletion skipped."
    fi
    echo
  else
    echo "No files to delete were found in directory $dir."
    echo
  fi
}

# Delete files in current directory
echo "Processing current directory..."
delete_files "."

# Check if test directory exists and process it
if [ -d "./test" ]; then
  echo "Processing ./test/ directory..."
  delete_files "./test"
else
  echo "Directory ./test/ does not exist."
fi

# Check if data directory exists and process it
if [ -d "./data" ]; then
  echo "Processing ./data/ directory..."
  delete_files "./data"
else
  echo "Directory ./data/ does not exist."
fi

# Check if data/bunny directory exists and process it
if [ -d "./data/bunny" ]; then
  echo "Processing ./data/bunny/ directory..."
  delete_files "./data/bunny"
else
  echo "Directory ./data/bunny/ does not exist."
fi
