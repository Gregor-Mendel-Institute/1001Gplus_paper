#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <new_directory>"
  exit 1
fi

new_directory="$1"
directories=("01_data" "02_scripts" "03_figures")

# Create the new directory if it doesn't exist
mkdir -p "$new_directory"

# Create subdirectories and README.md files inside the new directory
for dir in "${directories[@]}"; do
  dir_path="$new_directory/$dir"
  mkdir -p "$dir_path"
  echo "# $dir" > "$dir_path/README.md"
done

# Create README.md in the new directory
echo "# Base Directory" > "$new_directory/README.md"

echo "Directories and README.md files created in $new_directory."


