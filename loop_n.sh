#!/bin/bash

# Check if an argument was provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <ordering_file> <name_prefix>"
    exit 1
fi

# The first command line argument is the path to the ordering file
ordering_file="$1"
name_prefix="$2"

# Loop from 2 to 27
for n in {2..30}
do
    echo "Running with the first $n fasta files from $ordering_file..."

    # Read the first n lines of the ordering file, create an array of file paths
    readarray -t files < <(head -n "$n" "$ordering_file")

    # Construct the file paths
    file_paths=""
    for file in "${files[@]}"
    do
        file_paths+="$file "
    done

    # Trim the last extra space
    file_paths=$(echo $file_paths | sed 's/ $//')

    # Run your script with these paths
    time python3 jaccard_dit_wrapper.py $file_paths --output ${name_prefix}_$n --vg_output_dir dist_test/${name_prefix}__${n} --min_jaccard_init -w 512 -m 256

    echo "Completed run with $n fasta files."
done

