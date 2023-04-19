import os
import sys
import shutil

def copy_files(start, end, source_folder, destination_folder):
    if not os.path.exists(source_folder):
        print("Error: Source folder '{}' does not exist.".format(source_folder))
        sys.exit(1)

    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    for entry in os.listdir(source_folder):
        source_file_path = os.path.join(source_folder, entry)

        if "ToT.QE.PMTeff.txt" in entry:
            try:
                number = int(''.join(filter(str.isdigit, entry)))
                if number >= start and number <= end:
                    shutil.copy(source_file_path, destination_folder)
                    print("Copied '{}' to '{}'".format(entry, destination_folder))
            except ValueError:
                pass

    print("File copying completed.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: {} <start_number> <end_number> <source_folder>".format(sys.argv[0]))
        sys.exit(1)

    start = int(sys.argv[1])
    end = int(sys.argv[2])
    source_folder = sys.argv[3]
    destination_folder = "."

    copy_files(start, end, source_folder, destination_folder)

