#!/bin/bash

# Check if all required parameters are passed
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <base_directory> <folder_prefix> <start_replicate> <end_replicate> <band_width>"
    exit 1
fi

# Assign command line arguments to variables
BASE_DIR="$1"
FOLDER_PREFIX="$2"
START_REP="$3"
END_REP="$4"
BAND_WIDTH="$5"


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

	FASTA_FILES="${REP_DIR}/*.fasta"
        # Collect all rotated fasta files in the directory
        ROTATED_FASTAS_DIR="${REP_DIR}_run_${i}_w512_m256/rotated_fastas"
        ROTATED_FASTA_FILES="${ROTATED_FASTAS_DIR}/*.fasta.rotated"

        # Define output directories for jaccard_dit_wrapper.py outputs
        VG_OUTPUT_DIR="${REP_DIR}_run_${i}"
        OUTPUT_NAME="${REP_DIR}_jorder.txt"

        # Create output directory if it doesn't exist
        #mkdir -p "$VG_OUTPUT_DIR"

        # Command to run the jaccard_dit_wrapper.py script with time measurements
		echo "Running Jaccard Dit Wrapper:"
		/usr/bin/time -v python3 jaccard_dit_wrapper.py $FASTA_FILES --vg_output_dir $VG_OUTPUT_DIR --output $OUTPUT_NAME --orientation --min_jaccard_init -w $BAND_WIDTH -m 256 2>&1 | tee "$REP_DIR/vg_performance_stats_${i}.txt"
        
        # Concatenate all fasta files into one multifasta file
        MULTIFASTA="${ROTATED_FASTAS_DIR}/combined_${i}.mfa"
        cat ${ROTATED_FASTA_FILES} > ${MULTIFASTA}

        # Run MAFFT on the combined multifasta file
        MAFFT_OUTPUT="${REP_DIR}/combined_${i}_mafft_aligned.mfa"
        mafft --retree 2 --maxiterate 0 ${MULTIFASTA} > ${MAFFT_OUTPUT}

        # Run calculate_wentropy.py on the aligned fasta
        WENTROPY_OUTPUT="${REP_DIR}/output_wentropy_${i}.tsv"
        python3 calculate_wentropy.py ${MAFFT_OUTPUT} -o ${WENTROPY_OUTPUT}

        echo "Replicate $i processing complete."
    else
        echo "Directory $REP_DIR does not exist. Skipping..."
    fi
done

echo "All specified replicates processed."
