#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPTS_DIR="$SCRIPT_DIR/SimulatorOutputs"
export PROJECT="$SCRIPT_DIR/SceneCameraRay53/SceneCameraRay/SceneInterrogationin3DEnvironment.prj"
export SIM="$SCRIPT_DIR/Simulator_1.slx"
export MCITY="$SCRIPT_DIR/mcity03(1)(1).jpg"
export SCRIPT_DIR
export PYTHON_SCRIPTS_DIR


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


# Ensure the script exits on error
set -e

echo "RUNNING MATLAB"

gnome-terminal -- matlab -batch "mother_script"

MATLAB_PID=$!  # Captures the PID of the terminal process, not MATLAB


# Find the actual MATLAB process started by this terminal
MATLAB_PID=$(pgrep -f "matlab.*batch.*mother_script")

if [ -z "$MATLAB_PID" ]; then
    echo "Error: MATLAB process not found. Exiting."
    exit 1
fi

echo "MATLAB process PID: $MATLAB_PID"

while true; do
    source flag.sh
    echo "i am here"
    # Check if the environment variable is set
    if [ "$(printenv MATLAB_F)" = "true" ]; then
        echo "MATLAB script completed. Proceeding to close MATLAB."
        break
    fi
    sleep 1  # Check every 1 second
done

# Kill the MATLAB process if running
kill -9 $MATLAB_PID
echo "MATLAB process terminated."

if [ $? -eq 0 ]; then
    echo "MATLAB script completed successfully."
else
    echo "MATLAB script encountered an error."
    exit 1
fi

#run python scripts

echo "STARTING VDM"

cd $PYTHON_SCRIPTS_DIR

python3 new_inputgenerator.py
python3 VDM.py
python3 merege_csv.py

echo "RESULTS GENERATED, PLOTTING ACCELRATIONS."

python3 plotgraph.py

echo "EXITING SETUP, ROLLING BACK AND SAVING"

newfolderpath=$(create_next_dataset_folder)

cp "vdm_input.csv" "$newfolderpath"
cp "steering.csv" "$newfolderpath"
cp "output_vdm.csv" "$newfolderpath"

# Define the file to overwrite
file="$SCRIPT_DIR/flag.sh"

# Overwrite the file with new content
cat > $file <<EOF
export MATLAB_F="false"
EOF

# Make the file executable
chmod +x $file

echo "SIMULATION TERMINATED"




