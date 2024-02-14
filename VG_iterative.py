import os
import sys
import tempfile
from datetime import datetime
import random
import re
import subprocess
import glob
import json
import shutil
import argparse  # Import the argparse module

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process fasta file and generate graph.')
    parser.add_argument('fasta_filename', type=str, help='Fasta filename')
    parser.add_argument('--output_graph_name', type=str, default='output_graph', help='Output graph name (default: output_graph)')
    parser.add_argument('--output_dir', type=str, default=os.path.join(os.path.realpath(os.getcwd()), 'output_dir'), help='Output directory (default: ./output_dir)')
    return parser.parse_args()

def main():
    args = parse_args()  # Parse arguments

    # Convert relative paths to absolute paths
    args.fasta_filename = os.path.abspath(args.fasta_filename)
    args.output_dir = os.path.abspath(args.output_dir)

    print(datetime.now())
    print(os.path.realpath(os.getcwd()))

    print(f"Processing file: {args.fasta_filename}")
    print(f"Output directory: {args.output_dir}")

    # Create output directory if it does not exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Make temporary directory
    tmp_dir = tempfile.mkdtemp(prefix='vg-')
    print(f"Temporary directory: {tmp_dir}")

    os.chdir(tmp_dir)
    os.makedirs('fasta_files')
    os.makedirs('alignments')

    # Split fasta file into individual sequences and save to the fasta_files/ dir
    # Count number of sequences in fasta file
    if args.fasta_filename.endswith(".gz"):
        nr = subprocess.check_output(f"zcat {args.fasta_filename} | grep '>' -c", shell=True, text=True)
        subprocess.call(f"zcat {args.fasta_filename} | awk 'BEGIN {{n_seq=0;}} /^>/ {{if(n_seq%1==0){{file=sprintf(\"fasta_files/seq_%d.fa\",n_seq);}} print >> file; n_seq++; next;}} {{ print >> file; }}'", shell=True)
    else:
        nr = subprocess.check_output(f"grep '>' {args.fasta_filename} -c", shell=True, text=True)
        subprocess.call(f"cat {args.fasta_filename} | awk 'BEGIN {{n_seq=0;}} /^>/ {{if(n_seq%1==0){{file=sprintf(\"fasta_files/seq_%d.fa\",n_seq);}} print >> file; n_seq++; next;}} {{ print >> file; }}'", shell=True)

    nr = int(nr.strip())
    print("Number of sequences in total: ", nr)

    os.chdir("fasta_files")

    # Choose a random file/sequence
    random_file = random.choice(os.listdir())
    name = re.search(r'^>[a-zA-Z0-9_.]+', open(random_file).read()).group(0)[1:]
    print(f"Choosing file {random_file} ({name}) to build graph from")

    os.chdir("..")

    # Construct vg graph from the chosen file and circularize it
    subprocess.run(f"vg construct -M fasta_files/{random_file} > graph.vg", shell=True)
    subprocess.run(f"vg circularize -p {name} graph.vg > graph_circ.vg", shell=True)
    subprocess.run(f"vg stats -z graph_circ.vg", shell=True)

    os.remove("graph.vg")
    shutil.rmtree("fasta_files")

    # Further processing steps...

    # Move the final graph file to the output directory and perform final operations
    shutil.move(f"graph_circ.vg", os.path.join(args.output_dir, f"{args.output_graph_name}.vg"))

    # Perform clean-up and output conversion
    # For example, convert the file to ODGI and GFA formats, and remove the temporary directory

    print(datetime.now())

if __name__ == "__main__":
    main()
