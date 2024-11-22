import argparse
from collections import deque
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count
import os

def generate_kmer_deque(seq, k):
    kmers = deque()
    length = len(seq)
    wrapped_seq = seq + seq[:k-1]
    for i in range(length):
        kmers.append(wrapped_seq[i:i+k])
    return kmers, length

def process_segment(args):
    static_kmers, rotating_list, start, end, k = args
    best_score = 0
    best_rotation = start
    best_mask = None
    for i in range(start, end):
        rotated_list = rotating_list[i:] + rotating_list[:i]
        mask = [1 if static_kmers[j] == rotated_list[j] else 0 for j in range(len(static_kmers))]
        score = sum(mask)
        if score > best_score:
            best_score = score
            best_rotation = i
            best_mask = mask
    print(f"Process segment {start}-{end}: Best score {best_score}, Rotation {best_rotation}")
    return best_score, best_rotation, best_mask

def find_best_rotation_and_mask(static_kmers, rotating_kmers, num_processes, k):
    rotating_list = list(rotating_kmers)
    segment_length = len(rotating_kmers) // num_processes
    tasks = [(static_kmers, rotating_list, i * segment_length, (i + 1) * segment_length, k) for i in range(num_processes)]
    tasks[-1] = (static_kmers, rotating_list, tasks[-1][2], len(rotating_kmers), k)

    with Pool(processes=num_processes) as pool:
        results = pool.map(process_segment, tasks)

    best_score, best_rotation, best_mask = max(results, key=lambda x: x[0])
    print(f"Overall best rotation {best_rotation} with score {best_score}")
    return best_rotation, best_mask

def find_longest_contiguous_ones(mask):
    longest_run = 0
    current_run = 0
    start_index = 0
    for i, val in enumerate(mask):
        if val == 1:
            current_run += 1
            if current_run > longest_run:
                longest_run = current_run
                start_index = i - longest_run + 1
        else:
            current_run = 0
    return start_index, longest_run

def write_sequence(seq, file_path, rotation=None):
    if rotation is not None:
        seq = seq[rotation:] + seq[:rotation]
    record = SeqRecord(Seq(seq), id=os.path.basename(file_path), description="Processed sequence")
    SeqIO.write(record, file_path, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Align multiple circular DNA sequences to find shared regions.")
    parser.add_argument("file_list", type=str, help="Text file containing paths to FASTA files.")
    parser.add_argument("output_dir", type=str, help="Directory to save the processed FASTA files.")
    parser.add_argument("-k", "--kmer_size", type=int, default=7, help="Size of the k-mers to use.")
    parser.add_argument("-p", "--processes", type=int, default=cpu_count(), help="Number of processes to use.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load all sequences
    with open(args.file_list, 'r') as file:
        fasta_files = [line.strip() for line in file.readlines()]

    # Find the shortest sequence for reference
    sequences = [(SeqIO.read(f, "fasta").seq, f) for f in fasta_files]
    sequences.sort(key=lambda x: len(x[0]))
    ref_seq, ref_file = sequences[0]
    ref_kmers, _ = generate_kmer_deque(ref_seq, args.kmer_size)

    final_mask = [1] * len(ref_kmers)
    for seq, f in sequences[1:]:  # Skip the first as it's the reference
        rotating_kmers, _ = generate_kmer_deque(seq, args.kmer_size)
        rotation, mask = find_best_rotation_and_mask(list(ref_kmers), rotating_kmers, args.processes, args.kmer_size)
        print(f"Mask for {f}: {mask[:50]}...")# Print the first 50 elements of the mask for brevity       
        final_mask = [final_mask[i] & mask[i] for i in range(len(final_mask))]
        print(str(sum(final_mask)))

    start_index, longest_run = find_longest_contiguous_ones(final_mask)
    print(f"Longest contiguous region: Start={start_index}, Length={longest_run}")

    # Adjust sequences based on the longest shared region
    for seq, f in sequences:
        adjusted_seq = seq[start_index:] + seq[:start_index]
        write_sequence(adjusted_seq, os.path.join(args.output_dir, os.path.basename(f)))

if __name__ == "__main__":
    main()
