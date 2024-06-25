#!/bin/bash

# Check if all required parameters are passed
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <base_directory> <folder_prefix> <start_replicate> <end_replicate>"
    exit 1
fi

# Assign command line arguments to variables
BASE_DIR="$1"
FOLDER_PREFIX="$2"
START_REP="$3"
END_REP="$4"

# Ensure the base directory ends without a slash
BASE_DIR="${BASE_DIR%/}"

# Iterate over the range of specified replicates
for ((i=START_REP; i<=END_REP; i++))
do
    # Define the directory for the current replicate using the folder prefix
    REP_DIR="${BASE_DIR}/${FOLDER_PREFIX}_${i}"
    
    # Check if the directory exists
    if [ -d "$REP_DIR" ]; then
        echo "Processing replicate in $REP_DIR..."

        # Collect all fasta files in the directory
        FASTA_FILES="$REP_DIR/*.fasta"

        # Define output directories for jaccard_dit_wrapper.py outputs
        VG_OUTPUT_DIR="${REP_DIR}_run_${i}"
        OUTPUT_NAME="${REP_DIR}_${i}"

        # Create output directory if it doesn't exist
        mkdir -p "$VG_OUTPUT_DIR"

        # Command to run the python script with time measurements
        CMD="/usr/bin/time -v python3 jaccard_dit_wrapper.py $FASTA_FILES --vg_output_dir $VG_OUTPUT_DIR --output $OUTPUT_NAME --orientation --min_jaccard_init -w 512 -m 256"

        # Run the command and log performance metrics to a file
        $CMD 2>&1 | tee "$REP_DIR/usage_stats.txt"

        echo "Replicate $i processing complete."
    else
        echo "Directory $REP_DIR does not exist. Skipping..."
    fi
done

echo "All specified replicates processed."

