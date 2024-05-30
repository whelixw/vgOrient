#!/bin/bash

# Check if sufficient arguments were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <fasta_files_pattern> <output_dir_base>"
    exit 1
fi

# Assign arguments to variables
fasta_pattern=$1
output_dir_base=$2
result_file="$2.csv"

# Header for the CSV file
echo "Band_width,Max_node_size,RealTime,UserTime,SystemTime,Memory" > $result_file
#echo $1
#echo $2
# Arrays of values for width (-w) and max node size (-m)

widths=(256 512)
modes=(32 64 128 256 512)

# Find files matching the pattern
#readarray -t dataset_files < <(ls "$fasta_pattern")

# Check if files were found
#if [ ${#dataset_files[@]} -eq 0 ]; then
#    echo "No files matched the pattern $fasta_pattern."
#    exit 1
#fi

# Convert array to a space-separated string of quoted filenames
#dataset_arguments=""
#for file in "${dataset_files[@]}"; do
#    dataset_arguments+=\"$(realpath "$file")\"' '
#done

# Loop over each combination of width and mode
for w in "${widths[@]}"; do
    for m in "${modes[@]}"; do
        # Define dynamic output directory based on parameters
        vg_output_dir=${output_dir_base}

        echo "Running test with -w ${w} and -m ${m}"
        
        # Execute Python script with all files at once, quoting and escaping correctly
        /usr/bin/time -f "%e,%U,%S,%M" -o current_stats.txt --append \
	    python3 jaccard_dit_wrapper.py $fasta_pattern \
            --output $output_dir_base --vg_output_dir $vg_output_dir --min_jaccard_init -w $w -m $m

        # Extract runtime and memory usage from the output
        runtime_memory=$(tail -n 1 current_stats.txt)

        # Log the results into the CSV file
        echo "${w},${m},${runtime_memory}" >> $result_file

        echo "Test completed for -w ${w} and -m ${m}"
    done
done

echo "All tests completed. Results are stored in $result_file"

# Clean up temporary stats file
rm current_stats.txt
