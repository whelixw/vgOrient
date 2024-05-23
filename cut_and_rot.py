import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import sys
import os

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
    stripped_paths = {path: {node.rstrip('+-') for node in nodes} for path, nodes in paths.items()}
    shared_nodes = set.intersection(*stripped_paths.values())

    if not shared_nodes:
        raise ValueError("No shared nodes found in all paths")

    # Find the maximum length and all nodes that have this maximum length
    max_length = max(node_lengths[node] for node in shared_nodes)
    max_nodes = [node for node in shared_nodes if node_lengths[node] == max_length]

    # Check for multiple nodes with the same max size
    if len(max_nodes) > 1:
        print(f"Multiple nodes share the same maximum size of {max_length}: {', '.join(max_nodes)}")

    # We still return only one node and its sequence for further processing
    longest_shared_node = max_nodes[0]  # Arbitrarily choose the first one if there are multiple
    return longest_shared_node, node_sequences[longest_shared_node]

def rotate_sequence(seq, anchor):
    pos = seq.find(anchor)
    if pos == -1:
        raise ValueError(f"Anchor sequence '{anchor}' not found in the sequence")
    return seq[pos:] + seq[:pos]

def rotate_fastas(fasta_files, anchor_sequence, output_folder):
    for fasta_file in fasta_files:
        output_file = os.path.join(output_folder, os.path.basename(fasta_file) + ".rotated")
        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                rotated_seq = rotate_sequence(str(record.seq), anchor_sequence)
                record.seq = Seq(rotated_seq)
                SeqIO.write(record, output_handle, "fasta")

def main(gfa_file, fasta_files, output_folder):
    graph, node_sequences, paths = parse_vg_gfa(gfa_file)
    longest_node, anchor_sequence = find_longest_shared_node(paths, node_sequences)
    rotate_fastas(fasta_files, anchor_sequence, output_folder)
    print(f"Rotated FASTA files using the longest shared node {longest_node} ({anchor_sequence}) as the anchor.")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <gfa_file> <output_folder> <fasta_files...>")
        sys.exit(1)

    gfa_file = sys.argv[1]
    output_folder = sys.argv[2]
    fasta_files = sys.argv[3:]

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    main(gfa_file, fasta_files, output_folder)

