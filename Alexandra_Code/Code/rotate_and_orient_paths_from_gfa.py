import sys
import re
import textwrap

path_name_re = r"P\t([\w.]+)"
chars_to_remove = "\t*P"
node_id_to_seq = dict()
cutoff_node = str(sys.argv[3])

with open(sys.argv[2], "w") as fasta_file_out:
    with open(sys.argv[1], "r") as gfa_file_in:

        for line in gfa_file_in:

            # Save all node sequences to dict
            if line.startswith("S"):
                seq_line = line.split()
                node_id_to_seq[seq_line[1]] = seq_line[2]

            elif line.startswith("P"):

                path_line = line.strip()
                path_name = re.match(path_name_re, line).group(1)
                path_line = re.sub(path_name_re, "", path_line)

                for char in chars_to_remove:
                    path_line = path_line.replace(char, "")

                path_list = path_line.split(',')

                # List of tuples [(node, orientation), (node, orientation), ...]
                path_tuple = []
                for node_orient in path_list:
                    path_tuple.append((node_orient[:-1],node_orient[-1]))

                # Find index of cutoff node
                index_list = [index for index, item in enumerate(path_tuple) if item[0] == cutoff_node]
                if len(index_list) > 1:
                    print("Path", path_name, "goes through node", cutoff_node, "multiple times!")
                index = index_list[0]

                # Rotate path depending on node orientation
                if path_tuple[index][1] == '+':
                    rotated_path_tuple = path_tuple[index:] + path_tuple[:index]

                elif path_tuple[index][1] == '-':
                    rotated_path_tuple = path_tuple[index::-1] + path_tuple[-1:index:-1]

                rotated_seq = ""
                for node_orient in rotated_path_tuple:
                    rotated_seq += node_id_to_seq[node_orient[0]]


                # Write rotated sequences to FASTA file
                print(">" + path_name, file=fasta_file_out)
                print('\n'.join(textwrap.wrap(rotated_seq, 60)), file=fasta_file_out)