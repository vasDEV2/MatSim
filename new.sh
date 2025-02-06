#!/bin/bash


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPTS_DIR="$SCRIPT_DIR/SimulatorOutputs"

create_next_dataset_folder() {
    local prefix="DataSet"
    local max_number=0

    # Find all folders starting with "DataSet" and extract the number
    for folder in "$PYTHON_SCRIPTS_DIR"/${prefix}[0-9]*; do
        if [ -d "$folder" ]; then
            # Extract the numeric suffix
            folder_name=$(basename "$folder")
            number="${folder_name#$prefix}" # Remove the prefix
            if [[ "$number" =~ ^[0-9]+$ ]]; then # Check if it's a number
                max_number=$((number > max_number ? number : max_number))
            fi
        fi
    done

    # Increment the max number and create a new folder
    new_number=$((max_number + 1))
    new_folder="$PYTHON_SCRIPTS_DIR/${prefix}${new_number}"
    mkdir -p "$new_folder"
    echo "$new_folder"

}

create_next_dataset_folder
