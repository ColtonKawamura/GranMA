#!/bin/bash

# Define the target directories
target_dirs=("Width10" "Width20" "Width50" "Width100")

# Loop through each target directory and move matching files
for dir in "${target_dirs[@]}"; do
  # Create the directory if it doesn't exist
  mkdir -p "$dir"
  
  # Move files that match the pattern into the respective directory
  mv *"${dir}"* "$dir/"
done
