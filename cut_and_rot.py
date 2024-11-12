#!/usr/bin/env python
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import os
import argparse

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
                from_node, to_node = parts[1], parts[3]
                graph.add_edge(from_node, to_node)
            elif parts[0] == 'P':  # Path
                path_name = parts[1]
                path_nodes = parts[2].split(',')
                paths[path_name].extend(path_nodes)

    return graph, node_sequences, paths

def find_longest_shared_node(paths, node_sequences):
    stripped_paths = {path: [node.rstrip('+-') for node in nodes] for path, nodes in paths.items()}
    shared_nodes = set.intersection(*map(set, stripped_paths.values()))

    if not shared_nodes:
        raise ValueError("No shared nodes found in all paths")

    node_lengths = {node: len(node_sequences[node]) for node in shared_nodes}
    longest_shared_node = max(shared_nodes, key=lambda node: node_lengths[node])
    return longest_shared_node, node_sequences[longest_shared_node]


def rotate_sequence_to_node(seq, node_sequence):
    positions = [i for i in range(len(seq)) if seq.startswith(node_sequence, i)]
    if not positions:
        raise ValueError(f"Node sequence '{node_sequence}' not found in the sequence")
    if len(positions) > 1:
        raise ValueError(f"Multiple occurrences of node sequence '{node_sequence}' found in the sequence")

    node_start = positions[0]
    return seq[node_start:] + seq[:node_start]

def rotate_fastas(fasta_files, anchor_node_sequence, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for fasta_file in fasta_files:
        fasta_basename = os.path.basename(fasta_file)[:-6]  # Assume extension removal
        output_file = os.path.join(output_dir, f"{fasta_basename}.fasta.rotated")

        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                rotated_seq = rotate_sequence_to_node(str(record.seq), anchor_node_sequence)
                record.seq = Seq(rotated_seq)
                SeqIO.write(record, output_handle, "fasta")

def main(gfa_file, fasta_files, output_dir):
    _, node_sequences, paths = parse_vg_gfa(gfa_file)
    longest_node, anchor_node_sequence = find_longest_shared_node(paths, node_sequences)

    rotate_fastas(fasta_files, anchor_node_sequence, output_dir)
    print(f"Rotated FASTA files using the longest shared node {longest_node} as the anchor.")
    #print(f"FASTA files rotated using the anchor node sequence found in the longest shared node.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rotate FASTA files based on the longest shared node sequence.")
    parser.add_argument("gfa_file", help="The input GFA file.")
    parser.add_argument("fasta_files", nargs='+', help="List of input FASTA files.")
    parser.add_argument("--output_dir", default="rotated_fastas", help="Output directory for the rotated FASTA files.")
    args = parser.parse_args()

    main(args.gfa_file, args.fasta_files, args.output_dir)
