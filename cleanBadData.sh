#!/bin/bash

# Path to the file that contains the list of .mat files to move
input_file="bad_list.txt"

# Folder to move the files into
destination_folder="moved_files"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder"

# Read each line from the file
while IFS= read -r line; do
  # Remove the "Skipping file due to NaN values:" part if it exists
  file_path=$(echo "$line" | sed 's/^Skipping file due to NaN values: //')

  # Check if the resulting path is a .mat file starting with "2D"
  if [[ "$file_path" == *"2D"*".mat" ]]; then
    # Move the file if it exists
    if [ -f "$file_path" ]; then
      echo "Moving $file_path to $destination_folder"
      mv "$file_path" "$destination_folder/"
    else
      echo "File not found: $file_path"
    fi
  fi
done < "$input_file"
