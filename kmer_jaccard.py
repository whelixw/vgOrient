#!/usr/bin/env python
from Bio import SeqIO

from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import argparse
from pathlib import Path
import os

def reverse_fasta(fasta_file, record_id):
    reversed_flag = False
    with open(fasta_file, "r+") as file:
        records = list(SeqIO.parse(file, "fasta"))
        for record in records:
            if record.id == record_id:
                record.seq = record.seq.reverse_complement()  # Reverse the sequence
                reversed_flag = True
        file.seek(0)  # Reset file pointer to the beginning
        SeqIO.write(records, file, "fasta")
        file.truncate()  # Remove the rest of the original content after the new end
    return reversed_flag

def get_fasta_file_for_seq_id(fasta_files, seq_id):
    # Strip off '_reversed' if present to match against the sequence ID
    for fasta_file in fasta_files:
        file_stem = Path(fasta_file).stem
        while file_stem.endswith('_reversed'):
            file_stem = file_stem[:-9]  # Remove '_reversed' from the end
        if file_stem == seq_id:
            return fasta_file
    return None

def generate_kmers(seq, k):
    circular_seq_forward = seq + seq[:k - 1]
    forward_kmers = set(circular_seq_forward[i:i + k] for i in range(len(seq)))
    reverse_complement = str(Seq(seq).reverse_complement())
    circular_seq_reverse = reverse_complement + reverse_complement[:k - 1]
    reverse_kmers = set(circular_seq_reverse[i:i + k] for i in range(len(seq)))
    return forward_kmers, reverse_kmers


def generate_2mer_counts(seq):
    """Generate 2-mer counts for both forward and reverse (complement) orientations of a sequence."""
    k = 2
    forward_counts = {}
    reverse_counts = {}

    # Generate forward 2-mers
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in forward_counts:
            forward_counts[kmer] += 1
        else:
            forward_counts[kmer] = 1

    # Generate reverse complement of the sequence
    reverse_seq = str(Seq(seq).reverse_complement())
    # Generate reverse 2-mers
    for i in range(len(reverse_seq) - k + 1):
        kmer = reverse_seq[i:i + k]
        if kmer in reverse_counts:
            reverse_counts[kmer] += 1
        else:
            reverse_counts[kmer] = 1

    return forward_counts, reverse_counts

def calc_abs_diff(counts1, counts2):
    """Calculate the sum of absolute differences of 2-mer counts."""
    all_kmers = set(counts1.keys()).union(set(counts2.keys()))
    total_diff = sum(abs(counts1.get(k, 0) - counts2.get(k, 0)) for k in all_kmers)
    return total_diff


def calculate_jaccard_distance(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return 1 - (intersection / union)

def calculate_jaccard_distances(forward_kmers1, reverse_kmers1, forward_kmers2, reverse_kmers2):
    jaccard_cis = calculate_jaccard_distance(forward_kmers1, forward_kmers2)
    jaccard_trans = calculate_jaccard_distance(forward_kmers1, reverse_kmers2)
    return jaccard_cis, jaccard_trans

def order_files_by_distance(sequence_ids, orientations, cis_distances, trans_distances):
    # Calculate the effective distance for each sequence based on its orientation
    effective_distances = np.zeros(len(sequence_ids))
    for i, seq_id in enumerate(sequence_ids):
        for j in range(len(sequence_ids)):
            if i != j:
                if orientations[i] == orientations[j]:  # Both sequences have the same orientation
                    effective_distances[i] += cis_distances[i, j]
                else:  # Sequences have different orientations
                    effective_distances[i] += trans_distances[i, j]
    ordered_indices = np.argsort(effective_distances)
    return [sequence_ids[i] for i in ordered_indices]


def check_and_flip_orientations(sequence_ids, cis_distances, trans_distances):
    orientations = [False] * len(sequence_ids)  # Start with all sequences in original orientation
    total_improvement = True

    while total_improvement:
        total_improvement = False
        improved_this_round = True
        
        while improved_this_round:
            improved_this_round = False
            for seq_index in range(len(sequence_ids)):
                current_total_distance = sum(cis_distances[seq_index] if orientations[seq_index] else trans_distances[seq_index])
                potential_total_distance = sum(trans_distances[seq_index] if orientations[seq_index] else cis_distances[seq_index])
                
                if potential_total_distance < current_total_distance:
                    orientations[seq_index] = not orientations[seq_index]
                    improved_this_round = True
                    total_improvement = True
                    print(f"Flipping orientation of {sequence_ids[seq_index]} to {'reverse' if orientations[seq_index] else 'original'} based on improvement in Jaccard distance")

    return orientations

def simulated_annealing(sequence_ids, cis_distances, trans_distances, cooling_rate=0.99, initial_temp=100.0):
    import math
    import random

    orientations = [False] * len(sequence_ids)
    temp = initial_temp

    while temp > 0.1:
        seq_index = random.randint(0, len(sequence_ids) - 1)
        current_total_distance = sum(cis_distances[i][j] if orientations[i] == orientations[j] else trans_distances[i][j] for i in range(len(sequence_ids)) for j in range(len(sequence_ids)) if i != j)
        
        # Flip orientation temporarily
        orientations[seq_index] = not orientations[seq_index]
        new_total_distance = sum(cis_distances[i][j] if orientations[i] == orientations[j] else trans_distances[i][j] for i in range(len(sequence_ids)) for j in range(len(sequence_ids)) if i != j)

        if new_total_distance < current_total_distance or math.exp((current_total_distance - new_total_distance) / temp) > random.random():
            # Accept the new orientation
            #print(f"Flipping orientation of {sequence_ids[seq_index]} at temperature {temp}")
            pass
        else:
            # Revert back if not accepting the new configuration
            orientations[seq_index] = not orientations[seq_index]

        temp *= cooling_rate

    return orientations

def brute_force_optimal_orientation(sequence_ids, cis_distances, trans_distances):
    from itertools import product

    num_sequences = len(sequence_ids)
    best_score = float('inf')
    best_orientation = None

    # Iterate over all possible combinations of True (reversed) and False (original)
    for orientations in product([True, False], repeat=num_sequences):
        total_distance = sum(cis_distances[i][j] if orientations[i] == orientations[j] else trans_distances[i][j]
                             for i in range(num_sequences) for j in range(num_sequences) if i != j)

        if total_distance < best_score:
            best_score = total_distance
            best_orientation = orientations

    print(f"Optimal orientation found with score {best_score}")
    return list(best_orientation)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cluster sequences based on Jaccard distance of k-mers and optionally adjust orientations.')
    parser.add_argument('fasta_files', nargs='+', help='Input FASTA files')
    parser.add_argument('-k', '--kmer_size', type=int, default=11, help='K-mer size (default: 11)')
    parser.add_argument('-mcj', '--min_jaccard_init', action='store_true', help='Select initial sequence based on minimum cumulative Jaccard distance')
    parser.add_argument('-o', '--output_file', type=str, default='kmer_output.txt', help='Output file path (default: "kmer_output.txt")')
    parser.add_argument('-or','--orientation', action='store_true', help='Optimize sequence orientations based on Jaccard distances')
    parser.add_argument('-m', '--method', choices=['jaccard', 'abs_diff'], default='jaccard', help='Method to calculate distances: "jaccard" for Jaccard distance of k-mers or "abs_diff" for sum of absolute differences of 2-mers')
    #parser.add_argument('--log_dir', help='Directory to save the log file indicating reversed sequences. If not provided, uses the parent directory of the first FASTA file.')
    args = parser.parse_args()

    k = args.kmer_size
    fasta_files = args.fasta_files
    #log_dir = args.log_dir or str(Path(fasta_files[0]).parent)
    #log_file_path = os.path.join(log_dir, "reversed_sequences.log")
    kmers_per_sequence = {}

    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            kmers_per_sequence[record.id] = generate_kmers(str(record.seq), k)

    sequence_ids = list(kmers_per_sequence.keys())
    cis_distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    trans_distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    if args.method == 'abs_diff':
    # Calculate distances using absolute differences of 2-mer counts
        for i, j in combinations(range(len(sequence_ids)), 2):
            forward_counts_i, reverse_counts_i = generate_2mer_counts(str(kmers_per_sequence[sequence_ids[i]][0]))
            forward_counts_j, reverse_counts_j = generate_2mer_counts(str(kmers_per_sequence[sequence_ids[j]][0]))
        
            # Calculate distances for both forward-forward and forward-reverse (trans) scenarios
            cis_dist = calc_abs_diff(forward_counts_i, forward_counts_j)
            trans_dist = calc_abs_diff(forward_counts_i, reverse_counts_j)
        
            cis_distances[i, j] = cis_distances[j, i] = cis_dist
            trans_distances[i, j] = trans_distances[j, i] = trans_dist
    else:
        for i, j in combinations(range(len(sequence_ids)), 2):
            forward_kmers_i, reverse_kmers_i = kmers_per_sequence[sequence_ids[i]]
            forward_kmers_j, reverse_kmers_j = kmers_per_sequence[sequence_ids[j]]
            cis_dist, trans_dist = calculate_jaccard_distances(forward_kmers_i, reverse_kmers_i, forward_kmers_j, reverse_kmers_j)
            cis_distances[i, j] = cis_distances[j, i] = cis_dist
            trans_distances[i, j] = trans_distances[j, i] = trans_dist

    # Initialize the dictionary to store paths to reversed files
    #reversed_files = {}
    #orientations = [False] * len(sequence_ids)  # False indicates original orientation
    #if args.orientation:
    ##CHANGE ALGORITHMS HERE
    #orientations = check_and_flip_orientations(sequence_ids, cis_distances, trans_distances)
    orientations = simulated_annealing(sequence_ids, cis_distances, trans_distances)
    reversed_files = {}  # Initialize the dictionary to store paths to reversed files if needed
    # Write out the new FASTA files for sequences that need to be reversed
    for i, orientation in enumerate(orientations):
        seq_id = sequence_ids[i]
        fasta_file = get_fasta_file_for_seq_id(fasta_files, seq_id)
        if orientation and args.orientation:  # If the sequence needs to be reversed or restored
            reversed_file = reverse_fasta(fasta_file, seq_id)
            #print(f"{seq_id} reversed")
            reversed_files[seq_id] = fasta_file
            #print(f"{action} sequence {seq_id} saved to {reversed_file}")
        else:
            reversed_files[seq_id] = fasta_file

    #else:
    #    orientations = [False] * len(sequence_ids)  # All sequences remain in original orientation
    #    reversed_files = {seq_id: get_fasta_file_for_seq_id(fasta_files, seq_id) for seq_id in sequence_ids}
    
    # Decide the initial sequence based on cumulative Jaccard distances if the relevant argument is provided
if args.min_jaccard_init:
    total_jaccard_distances = np.sum(cis_distances, axis=1) + np.sum(trans_distances, axis=1)
    starting_index = np.argmin(total_jaccard_distances)
    ordered_sequence_ids = [sequence_ids[starting_index]]
    remaining_indices = list(range(len(sequence_ids)))
    remaining_indices.remove(starting_index)

    while remaining_indices:
        last_index = sequence_ids.index(ordered_sequence_ids[-1])  # This gets the integer index of the last ordered ID
        next_index = min(remaining_indices, key=lambda x: cis_distances[last_index, x] if orientations[x] == orientations[last_index] else trans_distances[last_index, x])
        ordered_sequence_ids.append(sequence_ids[next_index])
        remaining_indices.remove(next_index)       
else:
    ordered_sequence_ids = order_files_by_distance(sequence_ids, orientations, cis_distances, trans_distances)

# Use the specified output file name for writing results
output_file_path = args.output_file
with open(output_file_path, "w") as output_file:
    for seq_id in ordered_sequence_ids:
        fasta_path = reversed_files[seq_id]
        #print(fasta_path
        absolute_path = Path(fasta_path).absolute()  # Converts to absolute path
        output_file.write(f"{absolute_path}\n")
    print("Process completed. Check 'kmer_output.txt' for the ordered FASTA file paths.")

