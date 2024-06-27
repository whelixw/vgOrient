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
        FASTA_FILES="${REP_DIR}/*.fasta"

        # Concatenate all fasta files into one multifasta file
        MULTIFASTA="${REP_DIR}/combined_${i}.mfa"
        cat ${FASTA_FILES} > ${MULTIFASTA}

        # Define output file for MARS
        MARS_OUTPUT="${REP_DIR}/mars_output_${i}.mfa"

        # Time and run MARS on the combined multifasta file, output performance stats to a tsv
        /usr/bin/time -v mars -a DNA -i ${MULTIFASTA} -o ${MARS_OUTPUT} -m 0 -T 10 2>&1 | tee "${REP_DIR}/mars_performance_stats_${i}.txt"

        # Run MAFFT on the MARS output, time and log performance
        MAFFT_OUTPUT="${REP_DIR}/mars_mafft_aligned_${i}.mfa"
        mafft --retree 2 --maxiterate 0 ${MARS_OUTPUT} > ${MAFFT_OUTPUT} 

        # Run calculate_wentropy.py on the aligned fasta, time and log performance
        WENTROPY_OUTPUT="${REP_DIR}/output_wentropy_${i}.tsv"
        /usr/bin/time -v python3 calculate_wentropy.py ${MAFFT_OUTPUT} -o ${WENTROPY_OUTPUT} 2>&1 | tee "${REP_DIR}/wentropy_performance_${i}.tsv"

        echo "Replicate $i processing complete."
    else
        echo "Directory $REP_DIR does not exist. Skipping..."
    fi
done

echo "All specified replicates processed."
