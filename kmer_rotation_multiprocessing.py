import argparse
from collections import deque
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count
import os

def generate_kmer_deque(seq, k):
    """ Generate a deque of k-mers for the given sequence, treating it as circular. """
    kmers = deque()
    length = len(seq)
    wrapped_seq = seq + seq[:k-1]
    for i in range(length):
        kmers.append(wrapped_seq[i:i+k])
    return kmers

def process_segment(args):
    """ Process a segment of rotations and find the best score and corresponding rotation. """
    static_kmers, rotating_list, start, end = args
    best_score = 0
    best_rotation = start
    length_static = len(static_kmers)
    for i in range(start, end):
        rotated_list = rotating_list[i:] + rotating_list[:i]
        score = sum(1 for j in range(length_static) if static_kmers[j] == rotated_list[j])
        if score > best_score:
            best_score = score
            best_rotation = i
    return best_score, best_rotation

def find_best_rotation(static_kmers, rotating_kmers, num_processes):
    rotating_list = list(rotating_kmers)
    segment_length = len(rotating_kmers) // num_processes
    tasks = [(static_kmers, rotating_list, i * segment_length, (i + 1) * segment_length) for i in range(num_processes)]
    tasks[-1] = (static_kmers, rotating_list, tasks[-1][2], len(rotating_kmers))  # Handle remainder

    with Pool(processes=num_processes) as pool:
        results = pool.map(process_segment, tasks)

    best_score, best_rotation = max(results, key=lambda x: x[0])
    return best_rotation, best_score

def write_sequence(seq, file_path, rotation=None):
    """ Write the sequence, optionally rotated, to a specified file. """
    if rotation is not None:
        seq = seq[rotation:] + seq[:rotation]
    record = SeqRecord(Seq(seq), id=os.path.basename(file_path), description="Processed sequence")
    SeqIO.write(record, file_path, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Align multiple circular DNA sequences to a reference.")
    parser.add_argument("file_list", type=str, help="Text file containing paths to FASTA files.")
    parser.add_argument("output_dir", type=str, help="Directory to save the processed FASTA files.")
    parser.add_argument("-k", "--kmer_size", type=int, default=7, help="Size of the k-mers to use.")
    parser.add_argument("-p", "--processes", type=int, default=cpu_count(), help="Number of processes to use.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.file_list, 'r') as file:
        fasta_files = [line.strip() for line in file.readlines()]

    # Load the reference sequence
    ref_seq = SeqIO.read(fasta_files[0], "fasta").seq
    static_kmers = generate_kmer_deque(ref_seq, args.kmer_size)
    write_sequence(ref_seq, os.path.join(args.output_dir, os.path.basename(fasta_files[0])))  # Write reference as is

    # Process each query sequence
    for fasta_file in fasta_files[1:]:
        query_seq = SeqIO.read(fasta_file, "fasta").seq
        # Determine which sequence to rotate based on length
        if len(ref_seq) >= len(query_seq):
            rotating_kmers = generate_kmer_deque(ref_seq, args.kmer_size)
            static_kmers = generate_kmer_deque(query_seq, args.kmer_size)
            rotation, score = find_best_rotation(list(static_kmers), rotating_kmers, num_processes=args.processes)
            rotation = (len(ref_seq) - rotation) % len(ref_seq)  # Adjust rotation index
        else:
            static_kmers = generate_kmer_deque(ref_seq, args.kmer_size)
            rotating_kmers = generate_kmer_deque(query_seq, args.kmer_size)
            rotation, score = find_best_rotation(list(static_kmers), rotating_kmers, num_processes=args.processes)

        output_path = os.path.join(args.output_dir, os.path.basename(fasta_file))
        write_sequence(query_seq, output_path, rotation=rotation)

if __name__ == "__main__":
    main()
