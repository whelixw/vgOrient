import os
import sys
import argparse
import tempfile
import subprocess
from datetime import datetime
from pathlib import Path
import shutil

def setup_temporary_directory(prefix='vg-'):
    """Create and return a temporary directory."""
    tmp_dir = tempfile.mkdtemp(prefix=prefix)
    print(f"Temporary directory: {tmp_dir}")
    return Path(tmp_dir)

def prepare_output_directory(base_name, w, m, append_wm=True, specified_output_dir=None):
    """Prepare the output directory based on the base name and optional specified directory."""
    print(base_name, w, m, append_wm, specified_output_dir)
#    if append_wm:
#        specified_output_dir = f"{specified_output_dir}_w{w}_m{m}"
    output_dir = Path(specified_output_dir if specified_output_dir else f"{dir_name}_output").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir



def move_directory(src, dest):
    src = str(src)
    dest = str(dest)
    if os.path.exists(dest):
        # If the destination exists and is not empty, handle it appropriately:
        # For example, remove it first (if you want to replace it), merge contents, etc.
        # Here, we just print a message, but you might want to add actual handling logic.
        print(f"Destination {dest} already exists. Handling existing directory.")
        # shutil.rmtree(dest)  # Uncomment to remove the existing directory (be cautious)
    # Now perform the move
    shutil.move(src, dest)

def subprocess_command(cmd, output_file=None):
    """Run a subprocess command with optional output redirection."""
    if output_file:
        with open(output_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
    else:
        subprocess.run(cmd, check=True)

def read_sequence_data(fasta_file):
    """Read sequence data from a FASTA file, skipping the header."""
    with open(fasta_file, 'r') as file:
        next(file)  # Skip the header line
        sequence = ''.join(line.strip() for line in file)
    return sequence

def augment_graph(tmp_dir, graph_circ, gam_file, threads):
    """Handles the augmentation and re-indexing of the graph."""
    new_graph = tmp_dir / "new_graph.vg"
    subprocess_command(["vg", "augment", graph_circ, gam_file, "-i", "-S", "-t", threads], output_file=new_graph)
    # Move new graph to replace the old circular graph
    os.replace(new_graph, graph_circ)
    # Re-index the graph
    subprocess_command(["vg", "index", "-x", graph_circ.with_suffix(".xg"), str(graph_circ)])
    subprocess_command(["vg", "prune", "-k", "48", graph_circ], output_file=str(graph_circ) + "_pruned.vg")
    subprocess_command(["vg", "index", "-g", graph_circ.with_suffix(".gcsa"), "-Z", "400", str(graph_circ) + "_pruned.vg"])

def main():
    last_time = datetime.now()
    print(last_time)
    print(Path.cwd().resolve())

    parser = argparse.ArgumentParser(description='Script to process a list of FASTA files and generate output graphs.')
    parser.add_argument('fasta_list', help='Text file containing the list of FASTA files to process.')
    parser.add_argument('-o', '--output_dir_name', help='The name of the output directory.', default='output_graphs')
    parser.add_argument('-w', '--band_width', type=int, default=512, help='Band width for mapping.')
    parser.add_argument('-m', '--max_node_length', type=int, default=512, help='Max node length for constructing graph')
    parser.add_argument('--append_wm', action='store_true', default=False, help='Append w and m values to the output directory name.')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads for map and augment')

    args = parser.parse_args()
    fasta_list_filename = args.fasta_list
    output_dir = prepare_output_directory(Path.cwd(), args.band_width, args.max_node_length, args.append_wm, args.output_dir_name)

    with open(fasta_list_filename, 'r') as file:
        fasta_files = [line.strip() for line in file.readlines() if line.strip()]

    if not fasta_files:
        print("No FASTA files provided.")
        sys.exit(1)

    tmp_dir = setup_temporary_directory()
    os.chdir(tmp_dir)
    Path('alignments').mkdir(exist_ok=True)

    initial_fasta = fasta_files[0]
    base_name = os.path.splitext(os.path.basename(initial_fasta))[0]
    graph_circ = Path(f"{base_name}_graph_circ.vg")

    initial_output = tmp_dir / f"{base_name}_initial_output.vg"
    subprocess_command(["vg", "construct", "-r", initial_fasta, "-m", str(args.max_node_length), "-t", str(args.threads)], output_file=initial_output)
    circular_output = tmp_dir / f"{base_name}_circularized.vg"
    subprocess_command(["vg", "circularize", "-p", base_name, str(initial_output)], output_file=circular_output)
    os.replace(circular_output, graph_circ)  # Replace the original graph with the circularized version

    # Index the initial circular graph
    subprocess_command(["vg", "stats", "-z", str(graph_circ)])
    subprocess_command(["vg", "index", "-x", str(graph_circ.with_suffix(".xg")), str(graph_circ)])

    subprocess_command(["vg", "prune", "-k", "48", str(graph_circ)], output_file=str(graph_circ) + "_pruned.vg")
    subprocess_command(["vg", "index", "-g", str(graph_circ.with_suffix(".gcsa")), "-Z", "400", str(graph_circ) + "_pruned.vg"])

    # Map, augment, and update stats for other FASTA files
    for fasta_file in fasta_files[1:]:
        file_base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        print(f"mapping: {file_base_name}")
        print(f"time: {datetime.now()}, diff: {datetime.now()-last_time}")
        sequence_str = read_sequence_data(fasta_file)
        gam_file = Path('alignments') / f"{file_base_name}.gam"
        subprocess_command(["vg", "map", "-s", sequence_str, "-V", file_base_name, "-g", str(graph_circ.with_suffix(".gcsa")), "-x", str(graph_circ.with_suffix(".xg")), "-w", str(args.band_width), "-t", str(args.threads)], output_file=gam_file)
        augment_graph(tmp_dir, graph_circ, gam_file, str(args.threads))

    # Convert final VG to GFA
    final_gfa_path = output_dir / f"{base_name}_graph.gfa"
    final_odgi_path = output_dir / f"{base_name}_graph.odgi"
    #threads should be max threads by default, so no need to add -t
    #subprocess_command(["vg", "convert", "-f", str(graph_circ), "-t", str(args.threads)], output_file=final_gfa_path)
    subprocess_command(["vg", "convert", "-f", str(graph_circ)], output_file=final_gfa_path)

    #subprocess_command(["vg_1.44.0", "convert", "-o", str(graph_circ)], output_file=final_odgi_path)

    print(f"Graph processing completed at {datetime.now()}")
    print(f"All output files are saved in {output_dir}")

    new_dir_name = "intermediate"
    new_path = os.path.join(output_dir, new_dir_name)
    move_directory(tmp_dir, new_path)
    print(f"Temporary files are saved in {new_path}")

if __name__ == "__main__":
    main()
