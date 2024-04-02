#todo: use arparse instead of sys.argv to make script in line with VG iterative
from bdsg.bdsg import ODGI
import subprocess
import sys
from pathlib import Path

def get_number_of_paths_per_node(x):
    node_paths_tuple.append((graph.get_id(x), graph.get_step_count(x)))
    return True

def get_path_support(graph, node_id):
    paths = []
    graph.for_each_step_on_handle(graph.get_handle(node_id), lambda x: paths.append(graph.get_path_name(graph.get_path_handle_of_step(x))) or True)
    return paths

if len(sys.argv) < 2:
    print("Usage: script.py <path_to_odgi_file> [path_to_vg_file]")
    sys.exit(1)

graph = ODGI()
odgi_file_path = Path(sys.argv[1])
graph.deserialize(str(odgi_file_path))
base_path = odgi_file_path.parent
base_name = odgi_file_path.stem

# Infer VG file path or use the provided second argument
if len(sys.argv) == 3:
    vg_file_path = Path(sys.argv[2])
else:
    vg_file_path = base_path / f"{base_name}.vg"

print(vg_file_path)
print(base_name)
cmd = f"vg paths -v {vg_file_path} -L | wc -l"
paths = subprocess.run(cmd, shell=True, capture_output=True, text=True)
number_of_paths = int(paths.stdout.strip())
print(f"The number of paths in the graphs is: {number_of_paths}")

# For each node get the number of paths
node_paths_tuple = []
graph.for_each_handle(lambda x: get_number_of_paths_per_node(x))

# Get the nodes through which all paths go through
nodes_to_consider = []
for t in node_paths_tuple:
    if t[1] == number_of_paths:
        nodes_to_consider.append(t[0])

# or find the next best node
while nodes_to_consider == []:
    number_of_paths -= 1
    for t in node_paths_tuple:
        if t[1] == number_of_paths:
            nodes_to_consider.append(t[0])

# Check that each path only traverses the node once
final_nodes = []
for node in nodes_to_consider:
    paths = get_path_support(graph, node)
    if len(set(paths)) == number_of_paths:
        final_nodes.append(node)

# Get the length of the final nodes
node_length_dict = dict()
for node in final_nodes:
    handle = graph.get_handle(node)
    node_length_dict[node] = graph.get_length(handle)




# Only print nodes with the maximum length
max_length = max(node_length_dict.values())
print("Nodes with maximum length of", max_length, "bp saved to", base_path / f"{base_name}_cutoff_nodes.txt")
#open the file to write the results
with open(base_path / f"{base_name}_cutoff_nodes.txt", "w") as out_file:
    for item in node_length_dict.items():
        if item[1] == max_length:
            #print(item[0])
            out_file.write(f"{item[0]}\n")
    

