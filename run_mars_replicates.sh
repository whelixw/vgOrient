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

        # Define output file for MARS and log file for time
        MARS_OUTPUT="${REP_DIR}/mars_output_${i}.mfa"
        TIME_LOG="${REP_DIR}/mars_time_log_${i}.tsv"
        CMD_OUTPUT="${REP_DIR}/mars_cmd_output_${i}.txt"

        # Time and run MARS on the combined multifasta file
        ( /usr/bin/time -v mars -a DNA -i ${MULTIFASTA} -o ${MARS_OUTPUT} -m 0 -T 10 2> $TIME_LOG ) > $CMD_OUTPUT

        # Run MAFFT on the MARS output
        MAFFT_OUTPUT="${REP_DIR}/mars_mafft_aligned_${i}.mfa"
        TIME_LOG_MAFFT="${REP_DIR}/mafft_time_log_${i}.tsv"
        CMD_OUTPUT_MAFFT="${REP_DIR}/mafft_cmd_output_${i}.txt"

        ( /usr/bin/time -v mafft --retree 2 --maxiterate 0 ${MARS_OUTPUT} 2> $TIME_LOG_MAFFT ) > $MAFFT_OUTPUT > $CMD_OUTPUT_MAFFT

        # Run calculate_wentropy.py on the aligned fasta
        WENTROPY_OUTPUT="${REP_DIR}/output_wentropy_${i}.tsv"
        TIME_LOG_WENTROPY="${REP_DIR}/wentropy_time_log_${i}.tsv"
        CMD_OUTPUT_WENTROPY="${REP_DIR}/wentropy_cmd_output_${i}.txt"

        ( /usr/bin/time -v python3 calculate_wentropy.py ${MAFFT_OUTPUT} -o ${WENTROPY_OUTPUT} 2> $TIME_LOG_WENTROPY ) > $CMD_OUTPUT_WENTROPY

        echo "Replicate $i processing complete."
    else
        echo "Directory $REP_DIR does not exist. Skipping..."
    fi
done

echo "All specified replicates processed."
