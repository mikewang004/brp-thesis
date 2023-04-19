#!/bin/bash

# Check for correct number of arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <start_number> <end_number> <source_folder>"
  exit 1
fi

start=$1
end=$2
source_folder=$3
destination_folder="."

# Check if source folder exists
if [ ! -d "$source_folder" ]; then
  echo "Error: Source folder '$source_folder' does not exist."
  exit 1
fi

# Create destination folder if it doesn't exist
if [ ! -d "$destination_folder" ]; then
  mkdir "$destination_folder"
fi

# Loop through files in source folder
for file in "$source_folder"/*; do
  # Extract the filename
  filename=$(basename "$file")

  # Check if filename contains the integer range and the string "ToT.QE.PMTeff.txt"
  if [[ $filename == *ToT.QE.PMTeff.txt* ]]; then
    number=$(echo "$filename" | grep -oE '[0-9]+')
    if ((number >= start)) && ((number <= end)); then
      cp "$file" "$destination_folder"
      echo "Copied '$filename' to '$destination_folder'"
    fi
  fi
done

echo "File copying completed."

