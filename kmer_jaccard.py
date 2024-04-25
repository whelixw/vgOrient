from Bio import SeqIO
from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import argparse

def generate_kmers(seq, k):
    circular_seq_forward = seq + seq[:k - 1]
    forward_kmers = set(circular_seq_forward[i:i + k] for i in range(len(seq)))
    reverse_complement = str(Seq(seq).reverse_complement())
    circular_seq_reverse = reverse_complement + reverse_complement[:k - 1]
    reverse_kmers = set(circular_seq_reverse[i:i + k] for i in range(len(seq)))
    return forward_kmers, reverse_kmers

def calculate_jaccard_distance(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return 1 - (intersection / union)

def calculate_jaccard_distances(forward_kmers1, reverse_kmers1, forward_kmers2, reverse_kmers2):
    jaccard_cis = calculate_jaccard_distance(forward_kmers1, forward_kmers2)
    jaccard_trans = calculate_jaccard_distance(forward_kmers1, reverse_kmers2)
    return jaccard_cis, jaccard_trans

def flip_sequence(seq_index, orientations, cis_distances, trans_distances, current_distances):
    orientations[seq_index] = not orientations[seq_index]  # Flip the orientation
    for i in range(len(orientations)):
        if i != seq_index:
            if orientations[seq_index]:
                # If we're flipping to reverse, take the trans distance and put it in the current_distances
                current_distances[seq_index][i] = current_distances[i][seq_index] = trans_distances[seq_index][i] #double assignment skips lookup
            else:
                # If we're flipping back to original, take the cis distance and put it in the current_distances
                current_distances[seq_index][i] = current_distances[i][seq_index] = cis_distances[seq_index][i]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Determine sequence orientations based on Jaccard distances.')
    parser.add_argument('fasta_files', nargs='+', help='Input FASTA files')
    parser.add_argument('-k', '--kmer_size', type=int, default=11, help='K-mer size (default: 11)')
    args = parser.parse_args()

    k = args.kmer_size
    fasta_files = args.fasta_files
    kmers_per_sequence = {}

    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            kmers_per_sequence[record.id] = generate_kmers(str(record.seq), k)

    sequence_ids = list(kmers_per_sequence.keys())
    cis_distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    trans_distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    current_distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    orientations = [False] * len(sequence_ids)  # Start with all sequences in original orientation

    for i, j in combinations(range(len(sequence_ids)), 2):
        forward_kmers_i, reverse_kmers_i = kmers_per_sequence[sequence_ids[i]]
        forward_kmers_j, reverse_kmers_j = kmers_per_sequence[sequence_ids[j]]
        cis_dist, trans_dist = calculate_jaccard_distances(forward_kmers_i, reverse_kmers_i, forward_kmers_j, reverse_kmers_j)
        cis_distances[i, j] = cis_distances[j, i] = cis_dist
        trans_distances[i, j] = trans_distances[j, i] = trans_dist
        current_distances[i, j] = current_distances[j, i] = cis_dist  # Initially same as cis distances

    # Greedy algorithm for sequence orientation
    improved = True
    while improved:
        improved = False
        for seq_index in range(len(sequence_ids)):
            current_total_distance = sum(current_distances[seq_index])
            flip_sequence(seq_index, orientations, cis_distances, trans_distances, current_distances)
            new_total_distance = sum(current_distances[seq_index])

            if new_total_distance >= current_total_distance:
                # Flip back since there's no improvement
                flip_sequence(seq_index, orientations, cis_distances, trans_distances, current_distances)
            else:
                improved = True

    # Output the orientations
    for seq_id, orientation in zip(sequence_ids, orientations):
        print(f"{seq_id}: {'Reverse' if orientation else 'Original'}")
