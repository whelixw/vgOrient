from Bio import SeqIO

from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import argparse
from pathlib import Path
import os

def reverse_fasta(fasta_file, record_id):
    # Check if the file has already been reversed
    #print("Checking file:", fasta_file)
    if fasta_file.endswith('_reversed.fasta'):
        #print("test")
        # Determine the new file name by removing the "_reversed" suffix
        new_file = fasta_file.replace('_reversed.fasta', '.fasta')
        # Reverse the reversal (restore to original)
        action = "Restored"
    else:
        # Determine the new file name by adding the "_reversed" suffix
        new_file = fasta_file.replace('.fasta', '_reversed.fasta')
        action = "Reversed"
    
    with open(fasta_file, "r") as infile, open(new_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id == record_id:
                record.seq = record.seq.reverse_complement()  # Always reverse the sequence data
            SeqIO.write(record, outfile, "fasta")
    
    # Remove the old file
    os.remove(fasta_file)
    
    return new_file, action


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
    # Initial orientation: False indicates original orientation
    orientations = [False] * len(sequence_ids)
    proposed_flips = [False] * len(sequence_ids)  # Track proposed flips based on potential improvements

    # First, determine flips based on individual potential improvements
    for seq_index in range(len(sequence_ids)):
        current_total_distance = sum(cis_distances[seq_index] if orientations[seq_index] else trans_distances[seq_index])
        potential_total_distance = sum(trans_distances[seq_index] if orientations[seq_index] else cis_distances[seq_index])

        if potential_total_distance < current_total_distance:
            proposed_flips[seq_index] = True  # Mark for flipping based on individual checks

    # Count proposed flips to decide on the global strategy
    num_flips_proposed = sum(proposed_flips)
    num_no_flips = len(sequence_ids) - num_flips_proposed

    # If the number of proposed flips is less than those not flipped, flip those proposed
    # Otherwise, flip the others to minimize total reversals
    if num_flips_proposed < num_no_flips:
        # Apply flips only to those marked for flipping based on individual checks
        for i in range(len(orientations)):
            if proposed_flips[i]:
                orientations[i] = not orientations[i]  # Apply the flip
                print(f"Flipping orientation of {sequence_ids[i]} due to individual Jaccard improvement")
    else:
        # Flip all sequences not marked for flipping
        for i in range(len(orientations)):
            if not proposed_flips[i]:
                orientations[i] = not orientations[i]  # Apply the flip
                print(f"Flipping orientation of {sequence_ids[i]} to minimize total reversals")

    return orientations



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cluster sequences based on Jaccard distance of k-mers and optionally adjust orientations.')
    parser.add_argument('fasta_files', nargs='+', help='Input FASTA files')
    parser.add_argument('-k', '--kmer_size', type=int, default=11, help='K-mer size (default: 11)')
    parser.add_argument('-o', '--orientation', action='store_true', help='Optimize sequence orientations based on Jaccard distances')
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

    # Initialize the dictionary to store paths to reversed files
    #reversed_files = {}
    #orientations = [False] * len(sequence_ids)  # False indicates original orientation
    if args.orientation:
        orientations = check_and_flip_orientations(sequence_ids, cis_distances, trans_distances)
        reversed_files = {}  # Initialize the dictionary to store paths to reversed files if needed
    
        # Write out the new FASTA files for sequences that need to be reversed
        for i, orientation in enumerate(orientations):
            seq_id = sequence_ids[i]
            fasta_file = get_fasta_file_for_seq_id(fasta_files, seq_id)
            if orientation:  # If the sequence needs to be reversed or restored

                reversed_file, action = reverse_fasta(fasta_file, seq_id)
                
                reversed_files[seq_id] = reversed_file
                print(f"{action} sequence {seq_id} saved to {reversed_file}")
            else:
                reversed_files[seq_id] = fasta_file

    else:
        orientations = [False] * len(sequence_ids)  # All sequences remain in original orientation
        reversed_files = {seq_id: get_fasta_file_for_seq_id(fasta_files, seq_id) for seq_id in sequence_ids}
    
    # Output to kmer_output.txt
    with open("kmer_output.txt", "w") as output_file:
        ordered_sequence_ids = order_files_by_distance(sequence_ids, orientations, cis_distances, trans_distances)
        for seq_id in ordered_sequence_ids:
            fasta_path = reversed_files[seq_id]
            absolute_path = Path(fasta_path).absolute()  # Converts to absolute path
            output_file.write(f"{absolute_path}\n")
    
    print("Process completed. Check 'kmer_output.txt' for the ordered FASTA file paths.")

