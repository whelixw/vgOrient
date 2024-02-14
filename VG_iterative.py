import os
import sys
import argparse
import tempfile
from datetime import datetime
import random
import re
import subprocess
import glob
import json
import shutil
from pathlib import Path

print(datetime.now())
print(Path.cwd().resolve())

def setup_temporary_directory(prefix='vg-'):
    """Create and return a temporary directory."""
    tmp_dir = tempfile.mkdtemp(prefix=prefix)
    print(f"Temporary directory: {tmp_dir}")
    return Path(tmp_dir)

def prepare_output_directory(base_name, specified_output_dir=None):
    """Prepare the output directory based on the base name and optional specified directory."""
    output_dir = Path(specified_output_dir if specified_output_dir else f"{base_name}_output").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir

# Set up argument parser
parser = argparse.ArgumentParser(description='Script to process FASTA files and generate output graphs.')
parser.add_argument('fasta_filename', help='The FASTA file to process.')
parser.add_argument('-o', '--output_graph_name', help='The name of the output graph file. Defaults to the base name of the input file.', default=None)
parser.add_argument('-d', '--output_dir_name', help='The name of the output directory. Defaults to a directory named after the input file.', default=None)

# Parse arguments using pathlib for path arguments
args = parser.parse_args()
fasta_filename = Path(args.fasta_filename)
base_name = fasta_filename.stem
output_graph_name = args.output_graph_name if args.output_graph_name else base_name
output_dir = Path(args.output_dir_name if args.output_dir_name else f"{base_name}_output").resolve()


# Assign arguments to variables
fasta_filename = args.fasta_filename
base_name = os.path.splitext(os.path.basename(fasta_filename))[0]
output_graph_name = args.output_graph_name if args.output_graph_name else base_name
output_dir_name = args.output_dir_name if args.output_dir_name else os.path.join(os.path.realpath(os.getcwd()), base_name + "_output")



output_dir.mkdir(parents=True, exist_ok=True)




# Make temporary directory
tmp_dir = setup_temporary_directory()
output_dir = prepare_output_directory(base_name, args.output_dir_name)
os.chdir(tmp_dir)
Path('fasta_files').mkdir(exist_ok=True)
Path('alignments').mkdir(exist_ok=True)

# Split fasta file into individual sequences and save to the fasta_files/ dir
# Count number of sequences in fasta file
if fasta_filename.endswith(".gz"):
    nr = subprocess.check_output(f"zcat {fasta_filename} | grep '>' -c", shell=True, text=True)
    subprocess.call(f"zcat {fasta_filename} | awk 'BEGIN {{n_seq=0;}} /^>/ {{if(n_seq%1==0){{file=sprintf(\"fasta_files/seq_%d.fa\",n_seq);}} print >> file; n_seq++; next;}} {{ print >> file; }}'", shell=True)
else:
    nr = subprocess.check_output(f"grep '>' {fasta_filename} -c", shell=True, text=True)
    subprocess.call(f"cat {fasta_filename} | awk 'BEGIN {{n_seq=0;}} /^>/ {{if(n_seq%1==0){{file=sprintf(\"fasta_files/seq_%d.fa\",n_seq);}} print >> file; n_seq++; next;}} {{ print >> file; }}'", shell=True)

nr = int(nr.strip())
print("Number of sequences in total: ", nr)

os.chdir("fasta_files")

# Choose a random file/sequence
random_file = random.choice(os.listdir())
name = re.search(r'^>[a-zA-Z0-9_.]+', open(random_file).read()).group(0)[1:]
print("Choosing file {} ({}) to build graph from".format(random_file, name))

os.chdir("..")

# Construct vg graph from the chosen file and circularize it
subprocess.run(f"vg construct -M fasta_files/{random_file} > graph.vg", shell=True)
subprocess.run(f"vg circularize -p {name} graph.vg > graph_circ.vg", shell=True)
subprocess.run(f"vg stats -z graph_circ.vg", shell=True)

#os.remove("graph.vg")
#os.remove(f"fasta_files/{random_file}")
Path("graph.vg").unlink(missing_ok=True)
Path(f"fasta_files/{random_file}").unlink(missing_ok=True)
os.chdir("fasta_files")

for i in range(1, nr):
    identity = 0

    # Index graph
    subprocess.run(f"vg index -x {tmp_dir}/graph_circ.xg {tmp_dir}/graph_circ.vg", shell=True)  # XG index
    subprocess.run(f"vg prune -k 48 {tmp_dir}/graph_circ.vg > {tmp_dir}/graph_circ_pruned.vg", shell=True)
    subprocess.run(f"vg index -g {tmp_dir}/graph_circ.gcsa -Z 400 {tmp_dir}/graph_circ_pruned.vg", shell=True)  # GCSA index
    os.remove(f"{tmp_dir}/graph_circ_pruned.vg")

    for file in glob.glob("*.fa"):
        # Extract path name
        name = re.search(r'^>[a-zA-Z0-9_.]+', open(file).read()).group(0)[1:]

        # convert mitogenome to a string and align to graph
        with open(file) as f:
            next(f)
            mitogenome_str = f.read().replace('\n', '')
        subprocess.run(f"vg map -s {mitogenome_str} -V {name} -g {tmp_dir}/graph_circ.gcsa -x {tmp_dir}/graph_circ.xg > {tmp_dir}/alignments/{file}.gam", shell=True)

        # Run vg view command and capture output
        command = f"vg view -a {tmp_dir}/alignments/{file}.gam"
        output = subprocess.check_output(command, shell=True)

        # Extract identity using jq
        data = json.loads(output)
        try:
            identity_new = data["identity"] # when encountering a file with no identity it crashes
         # Set identity to 0.0 if it is "null"
        except:
            identity_new = 0.0

        # Save the file with the highest %identity
        if float(identity_new) >= identity:
            identity = float(identity_new)
            best_file = file
            best_name = name

        print(f"{file} ({name}) aligned with identity score {identity_new}")

    print(f"{best_file} has the highest %identity {identity}")

  # Augment graph with that sequence
    subprocess.run(f"vg augment -i -S {tmp_dir}/graph_circ.vg {tmp_dir}/alignments/{best_file}.gam  > {tmp_dir}/graph_aug.vg", shell=True)

    #subprocess.run(f"rm {tmp_dir}/fasta_files/{best_file}", shell=True)
    os.remove(os.path.join(tmp_dir, "fasta_files", best_file))
    subprocess.run(f"mv {tmp_dir}/graph_aug.vg {tmp_dir}/graph_circ.vg", shell=True)

    subprocess.run(f"vg stats -z {tmp_dir}/graph_circ.vg", shell=True)

    # Clear the alignments dir for the next round
    #could do rmtree followed by mkdir instead
    #for file in os.listdir(os.path.join(tmp_dir, "alignments")):
        #os.remove(os.path.join(tmp_dir, "alignments", file))
    for file in Path(tmp_dir).joinpath("alignments").iterdir():
        file.unlink(missing_ok=True)

# Move the graph file to the output directory
#shutil.move(f"{tmp_dir}/graph_circ.vg", f"{output_dir_name}/{output_graph_name}.vg")
#shutil.move(os.path.join(tmp_dir, "graph_circ.vg"), os.path.join(output_dir_name, f"{output_graph_name}.vg"))
#shutil.move(Path(tmp_dir).joinpath("graph_circ.vg"), output_dir.joinpath(f"{output_graph_name}.vg"))

# Change to the new directory
os.chdir(output_dir_name)

# Convert the file to ODGI and GFA formats
subprocess.run(f"vg_1.44.0 convert -o {output_graph_name}.vg > {output_graph_name}.odgi", shell=True)

subprocess.run(f"vg convert -f {output_graph_name}.vg > {output_graph_name}.gfa", shell=True)

# Remove the temporary directory
os.system(f"rm -rf {tmp_dir}")

# Print the current date and time
print(datetime.now())


