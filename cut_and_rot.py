import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import os
import sys

def parse_vg_gfa(gfa_file):
    graph = nx.DiGraph()
    node_sequences = {}
    paths = defaultdict(list)

    with open(gfa_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'S':  # Segment
                node_id = parts[1]
                sequence = parts[2]
                node_sequences[node_id] = sequence
                graph.add_node(node_id, sequence=sequence)
            elif parts[0] == 'L':  # Link
                from_node = parts[1]
                to_node = parts[3]
                graph.add_edge(from_node, to_node)
            elif parts[0] == 'P':  # Path
                path_name = parts[1]
                path_nodes = parts[2].split(',')
                paths[path_name].extend(path_nodes)
    
    return graph, node_sequences, paths

def find_longest_shared_node(paths, node_sequences):
    node_lengths = {node: len(seq) for node, seq in node_sequences.items()}
    # Strip orientation symbols from node identifiers in paths
    stripped_paths = {path: {node.rstrip('+-') for node in nodes} for path, nodes in paths.items()}
    shared_nodes = set.intersection(*stripped_paths.values())

    if not shared_nodes:
        raise ValueError("No shared nodes found in all paths")

    longest_shared_node = max(shared_nodes, key=lambda node: node_lengths[node])
    return longest_shared_node, node_sequences[longest_shared_node]

def rotate_sequence_to_node(seq, node_sequence):
    node_start = seq.find(node_sequence)
    if node_start == -1:
        raise ValueError(f"Node sequence '{node_sequence}' not found in the sequence")
    
    return seq[node_start:] + seq[:node_start]

def rotate_fastas(fasta_files, anchor_node_sequence, output_files):
    for fasta_file, output_file in zip(fasta_files, output_files):
        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                rotated_seq = rotate_sequence_to_node(str(record.seq), anchor_node_sequence)
                record.seq = Seq(rotated_seq)
                SeqIO.write(record, output_handle, "fasta")

def main(gfa_file, fasta_files, output_dir):
    graph, node_sequences, paths = parse_vg_gfa(gfa_file)
    longest_node, anchor_node_sequence = find_longest_shared_node(paths, node_sequences)
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    output_files = [os.path.join(output_dir, os.path.basename(fasta_file) + ".rotated") for fasta_file in fasta_files]

    rotate_fastas(fasta_files, anchor_node_sequence, output_files)
    print(f"Rotated FASTA files using the longest shared node {longest_node} ({anchor_node_sequence}) as the anchor.")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <gfa_file> <output_directory> <fasta_files...>")
        sys.exit(1)

    gfa_file = sys.argv[1]
    output_dir = sys.argv[2]
    fasta_files = sys.argv[3:]

    main(gfa_file, fasta_files, output_dir)
