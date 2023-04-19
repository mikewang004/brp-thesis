#!/bin/bash

# Usage: ./move_files.sh <start_number> <end_number> <source_folder>

start=14399
end=14422
source_folder=../../../../../sps/km3net/repo/data/calibration/KM3NeT_00000133/v8.0_PMTeff_new/calibration/
destination_folder=$(dirname "$0")/

if [ ! -d "$destination_folder" ]; then
  echo "Error: Destination folder does not exist."
  exit 1
fi

if [ ! -d "$source_folder" ]; then
  echo "Error: Source folder does not exist."
  exit 1
fi

for file in "$source_folder"/*ToT.QE.PMTeff.txt; do
  filename=$(basename "$file")
  number=$(echo "$filename" | grep -oE '[0-9]+')
  if [ -n "$number" ] && (( number >= start && number <= end )); then
    mv "$file" "$destination_folder"
    echo "Moved $filename to $destination_folder"
  fi
done

echo "File movement completed."
