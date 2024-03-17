#!/usr/bin/env python3
import argparse
class GFA:
    def __init__(self, file_path):
        self.dict_of_node_length = {}
        self.paths = {}
        self.load_gfa(file_path)
        self.name = file_path

    def load_gfa(self, file_path):
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("S"):
                    fields = line.strip().split("\t")
                    node_id = fields[1]
                    sequence = fields[2]
                    self.dict_of_node_length[node_id] = len(sequence)
                elif line.startswith("P"):
                    fields = line.strip().split("\t")
                    path_name = fields[1]
                    node_ids = fields[2].split(",")
                    node_ids = [node_id[:-1] for node_id in node_ids]  # Removes the last character (+ or -)
                    self.paths[path_name] = node_ids

    def find_node_positions(self, target_node, target_path=None):
        results = []
        target_node = str(target_node)
        for path_name, node_ids in self.paths.items():
            if target_path and path_name != target_path:
                continue
            current_length = 0
            for node_id in node_ids:
                node_length = self.dict_of_node_length[node_id]
                if node_id == target_node:
                    start_pos = current_length
                    end_pos = current_length + node_length
                    results.append((path_name, node_id, start_pos, end_pos))
                current_length += node_length
        return results

    def find_nodes_in_range(self, start_idx, end_idx, target_path):
        if target_path not in self.paths:
            return []
        nodes_in_range = []
        current_length = 0
        for node_id in self.paths[target_path]:
            node_length = self.dict_of_node_length[node_id]
            node_end = current_length + node_length
            if node_end > start_idx and current_length < end_idx:
                nodes_in_range.append(node_id)
            current_length += node_length
        return nodes_in_range


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GFA files.")
    parser.add_argument("gfa_paths", nargs=2, help="Paths to two GFA files.")
    parser.add_argument("--node", "-n", help="Target node ID to find in all paths.")
    parser.add_argument("--path", "-p", help="Target path name for finding nodes in a range.")
    parser.add_argument("--start", "-s", type=int, help="Start index for range query.")
    parser.add_argument("--end", "-e", type=int, help="End index for range query.")

    args = parser.parse_args()

    # Load GFA files
    gfa1 = GFA(args.gfa_paths[0])
    gfa2 = GFA(args.gfa_paths[1])

    # Node positions query
    if args.node:
        for gfa in [gfa1, gfa2]:
            positions = gfa.find_node_positions(args.node, args.path)
            for position in positions:
                print(f"File: {gfa.name}\nPath: {position[0]}, Node: {position[1]}, Start: {position[2]}, End: {position[3]}")

    # Nodes in range query
    if args.path and args.start is not None and args.end is not None:
        for gfa in [gfa1, gfa2]:
            nodes = gfa.find_nodes_in_range(args.start, args.end, args.path)
            print(f"File: {gfa.name}\nNodes in range {args.start}-{args.end} for path {args.path} in GFA: {nodes}")