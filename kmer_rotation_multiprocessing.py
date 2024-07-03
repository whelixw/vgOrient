import argparse
from collections import deque
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

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
    for i in range(start, end):
        rotated_list = rotating_list[i:] + rotating_list[:i]  # Rotate the k-mers
        score = sum(1 for j, kmer in enumerate(static_kmers) if kmer == rotated_list[j])
        if score > best_score:
            best_score = score
            best_rotation = i
    return best_score, best_rotation

def find_best_rotation(static_kmers, rotating_kmers, num_processes):
    """ Determine the best rotation using multiprocessing. """
    rotating_list = list(rotating_kmers)
    segment_length = len(rotating_kmers) // num_processes
    tasks = [(static_kmers, rotating_list, i * segment_length, (i + 1) * segment_length) for i in range(num_processes)]
    
    # Handle any remaining rotations in the last segment
    tasks[-1] = (static_kmers, rotating_list, tasks[-1][2], len(rotating_kmers))
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_segment, tasks)
    
    # Find the best result across all segments
    best_score, best_rotation = max(results, key=lambda x: x[0])
    return best_rotation, best_score

def main():
    parser = argparse.ArgumentParser(description="Align two circular DNA sequences using multiprocessing.")
    parser.add_argument("fasta1", type=str, help="FASTA file containing the first sequence.")
    parser.add_argument("fasta2", type=str, help="FASTA file containing the second sequence.")
    parser.add_argument("-k", "--kmer_size", type=int, default=7, help="Size of the k-mers to use.")
    parser.add_argument("-p", "--processes", type=int, default=cpu_count(), help="Number of processes to use.")

    args = parser.parse_args()

    seq1 = SeqIO.read(args.fasta1, "fasta").seq
    seq2 = SeqIO.read(args.fasta2, "fasta").seq
    k = args.kmer_size
    processes = args.processes

    if len(seq1) >= len(seq2):
        static_kmers = generate_kmer_deque(seq2, k)
        rotating_kmers = generate_kmer_deque(seq1, k)
    else:
        static_kmers = generate_kmer_deque(seq1, k)
        rotating_kmers = generate_kmer_deque(seq2, k)

    rotation, score = find_best_rotation(list(static_kmers), rotating_kmers, num_processes=processes)
    print(f"Best rotation: {rotation} positions with an alignment score of {score}")

if __name__ == "__main__":
    main()
