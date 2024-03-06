import sys

def read_fasta(file_path):
    """Read a FASTA file and return a dictionary of sequences indexed by their headers."""
    sequences = {}
    header = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                header = line
                sequences[header] = []
            else:
                sequences[header].append(line)
    # Join lines to form single continuous strings for each sequence
    for header in sequences:
        sequences[header] = ''.join(sequences[header])
    return sequences

def write_fasta(sequences, output_path, wrap=80):
    """Write sequences to a FASTA file with specified character wrap."""
    with open(output_path, 'w') as output_file:
        for header, sequence in sequences.items():
            output_file.write(header + '\n')
            for i in range(0, len(sequence), wrap):
                output_file.write(sequence[i:i+wrap] + '\n')

def offset_sequence(sequence, n):
    """Offset a sequence by n characters, treating it as circular."""
    offset = n % len(sequence)  # Handle offsets larger than the sequence
    return sequence[offset:] + sequence[:offset]

def main(input_file, output_file, offset_n):
    sequences = read_fasta(input_file)
    offsetted_sequences = {header: offset_sequence(seq, int(offset_n)) for header, seq in sequences.items()}
    write_fasta(offsetted_sequences, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_fasta> <output_fasta> <offset>")
        sys.exit(1)
    input_file, output_file, offset = sys.argv[1], sys.argv[2], sys.argv[3]
    main(input_file, output_file, offset)

