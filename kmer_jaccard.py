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

def adjust_orientations(sequence_ids, kmers_per_sequence, cis_distances, trans_distances):
    orientations = {seq_id: 'original' for seq_id in sequence_ids}  # All sequences start in the original orientation
    improved = True

    while improved:
        improved = False
        for i in range(len(sequence_ids)):
            current_orientation = orientations[sequence_ids[i]]
            # Calculate current total distance with the current orientation
            current_total_distance = np.sum(cis_distances[i] if current_orientation == 'original' else trans_distances[i])
            
            # Calculate potential total distance with the flipped orientation
            potential_total_distance = np.sum(trans_distances[i] if current_orientation == 'original' else cis_distances[i])
            
            # Decide to flip if it reduces the total distance
            if potential_total_distance < current_total_distance:
                print(potential_total_distance)
                print(current_total_distance)
                print(cis_distances[i])
                print(trans_distances[i])
                orientations[sequence_ids[i]] = 'reverse' if current_orientation == 'original' else 'original'
                improved = True
                # Recalculate the cis and trans distances for this sequence against all others

    return orientations

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cluster sequences based on Jaccard distance of k-mers.')
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

    for i, j in combinations(range(len(sequence_ids)), 2):
        forward_kmers_i, reverse_kmers_i = kmers_per_sequence[sequence_ids[i]]
        forward_kmers_j, reverse_kmers_j = kmers_per_sequence[sequence_ids[j]]
        cis_dist, trans_dist = calculate_jaccard_distances(forward_kmers_i, reverse_kmers_i, forward_kmers_j, reverse_kmers_j)
        cis_distances[i, j] = cis_distances[j, i] = cis_dist
        trans_distances[i, j] = trans_distances[j, i] = trans_dist

    # Adjust orientations based on the Jaccard distances
    final_orientations = adjust_orientations(sequence_ids, kmers_per_sequence, cis_distances, trans_distances)
    print("Final Orientations:", final_orientations)
