class GFA:
    def __init__(self, file_path):
        self.dict_of_node_length = {}
        self.paths = {}
        self.load_gfa(file_path)

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


# Load two GFA files into objects
gfa1 = GFA("suina_2_MW534270.1.gfa")
gfa2 = GFA("suina_2_DQ409327.1.gfa")

# Example usage:

# Find the start and end positions of a target node in all paths for both GFAs
target_node = "3569"
for gfa in [gfa1, gfa2]:
    positions = gfa.find_node_positions(target_node)
    for position in positions:
        print(f"GFA: {position[0]}, Node: {position[1]}, Start: {position[2]}, End: {position[3]}")

# Given a range of indices and a path, retrieve the nodes covered by those indices
start_index = 14352
end_index = 14456
target_path = "DQ409327.1"
for gfa in [gfa1, gfa2]:
    nodes = gfa.find_nodes_in_range(start_index, end_index, target_path)
    print(f"Nodes in range {start_index}-{end_index} for path {target_path} in GFA: {nodes}")
