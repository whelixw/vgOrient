from Bio import SeqIO
from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import argparse


def generate_kmers(seq, k):
    """
    Generate k-mers for both forward and reverse (complement) orientations of a circular sequence.

    :param seq: The sequence (as a string).
    :param k: The k-mer size.
    :return: A tuple containing two sets of k-mers: (forward_kmers, reverse_kmers)
    """
    # Extend the sequence for circularity in the forward direction
    circular_seq_forward = seq + seq[:k - 1]
    forward_kmers = set(circular_seq_forward[i:i + k] for i in range(len(seq)))

    # Generate the reverse complement and extend for circularity in the reverse direction
    reverse_complement = str(Seq(seq).reverse_complement())
    circular_seq_reverse = reverse_complement + reverse_complement[:k - 1]
    reverse_kmers = set(circular_seq_reverse[i:i + k] for i in range(len(seq)))

    return forward_kmers, reverse_kmers


def calculate_jaccard_distance(set1, set2):
    """
    Calculate the Jaccard distance between two sets.

    :param set1: The first set of k-mers.
    :param set2: The second set of k-mers.
    :return: The Jaccard distance.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return 1 - (intersection / union)


def average_jaccard(debug, forward_kmers1, reverse_kmers1, forward_kmers2, reverse_kmers2):
    if debug:
        print("Forward 1 Sample: ", list(forward_kmers1)[:5])
        print("Reverse 1 Sample: ", list(reverse_kmers1)[:5])
        print("Forward 2 Sample: ", list(forward_kmers2)[:5])
        print("Reverse 2 Sample: ", list(reverse_kmers2)[:5])
    """
    Calculate the average Jaccard distance between two sequences, considering both forward and reverse directions.

    :param forward_kmers1: Forward k-mers of the first sequence.
    :param reverse_kmers1: Reverse k-mers of the first sequence.
    :param forward_kmers2: Forward k-mers of the second sequence.
    :param reverse_kmers2: Reverse k-mers of the second sequence.
    :return: The average Jaccard distance.
    """
    # Calculate Jaccard distances in both directions
    jaccard_cis = calculate_jaccard_distance(forward_kmers1, forward_kmers2)
    jaccard_trans = calculate_jaccard_distance(forward_kmers1, reverse_kmers2)

    # Return the average of the two distances
    if jaccard_cis == jaccard_trans:
        if debug:
            print("cis and trans Jaccard distances are equal")
        return jaccard_cis
    else:
        if debug:
            print("cis and trans Jaccard distances are not equal:", jaccard_cis, jaccard_trans)
    return (jaccard_cis + jaccard_trans) / 2





def calculate_distance_vector_and_min_index(distances, cluster_indices):
    """
    Calculate the distance vector for the given cluster indices and find the index of the minimum value.

    :param distances: A numpy array representing the distance matrix.
    :param cluster_indices: A list of indices representing the current cluster.
    :return: A numpy array representing the distance vector and the index of the minimum value in the original distances matrix.
    """
    # Initialize the distance vector with NaNs, which will be ignored in the mean calculation
    distance_vector = np.full(distances.shape[1], np.nan)

    for i in range(distances.shape[1]):
        # If the column index is in cluster_indices, skip it as we don't calculate distance to itself
        if i in cluster_indices:
            continue

        # Select rows corresponding to the cluster indices and the current column
        values = distances[cluster_indices, i]

        # Check if any value is zero, and if not, calculate the mean for the column
        if not np.any(values == 0):
            distance_vector[i] = np.mean(values)

    # Find the index of the minimum value in the distance vector, considering only non-NaN values
    if np.all(np.isnan(distance_vector)):  # Check if all values are NaN
        min_index = None  # No valid minimum if all values are NaN
    else:
        min_index = np.nanargmin(distance_vector)  # Find index of the minimum non-NaN value

    # Return both the cleaned-up distance vector (with NaNs removed) and the index of the minimum value
    return distance_vector[~np.isnan(distance_vector)], min_index


def find_closest_sequences(distances):
    """
    Find the two closest sequences based on the distance matrix.

    :param distances: A numpy array representing the distance matrix.
    :return: A tuple of indices representing the two closest sequences.
    """
    # Mask the diagonal to ignore self-distances
    np.fill_diagonal(distances, np.inf)

    # Find the indices of the minimum value in the distance matrix
    min_index = np.unravel_index(np.argmin(distances, axis=None), distances.shape)

    # Restore the diagonal back to 0
    np.fill_diagonal(distances, 0)

    return min_index







def iterative_clustering(distances):
    """
    Perform iterative clustering by starting with the closest two sequences and then
    adding sequences to the cluster based on the minimum distance to the cluster.

    :param distances: A numpy array representing the distance matrix.
    :return: The order in which sequences are added to the cluster.
    """
    # Step 1: Find the initial two closest sequences
    i, j = find_closest_sequences(distances)
    cluster_indices = [i, j]  # Initialize the cluster with these two sequences

    # Keep track of the order in which indices are added to the cluster
    clustering_order = [i, j]

    # Step 2: Iteratively add sequences to the cluster
    while len(cluster_indices) < distances.shape[0]:
        # Calculate the distance vector for the current cluster and find the next sequence to add
        _, min_index = calculate_distance_vector_and_min_index(distances, cluster_indices)

        # If min_index is None, it means all sequences have been considered; break the loop
        if min_index is None:
            break

        # Add the sequence with the minimum distance to the cluster_indices
        if min_index not in cluster_indices:
            cluster_indices.append(min_index)
            clustering_order.append(min_index)

    return clustering_order

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Cluster sequences based on Jaccard distance of k-mers.')
    parser.add_argument('fasta_files', nargs='+', help='Input FASTA files')
    parser.add_argument('-k', '--kmer_size', type=int, default=11, help='K-mer size (default: 11)')
    args = parser.parse_args()

    # Store the k-mer size and input files
    k = args.kmer_size
    fasta_files = args.fasta_files
    kmers_per_sequence = {}

    # Generate k-mers for each sequence
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            kmers_per_sequence[record.id] = generate_kmers(str(record.seq), k)

    # Your existing code for calculating distances and clustering
    # Replace the hardcoded fasta_files and k values with the parsed arguments
    sequence_ids = list(kmers_per_sequence.keys())
    distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    for i, j in combinations(range(len(sequence_ids)), 2):
        # Unpack the tuples of forward and reverse k-mers for each sequence
        forward_kmers_i, reverse_kmers_i = kmers_per_sequence[sequence_ids[i]]
        forward_kmers_j, reverse_kmers_j = kmers_per_sequence[sequence_ids[j]]
        # Calculate the average Jaccard distance
        dist = average_jaccard(False, forward_kmers_i, reverse_kmers_i, forward_kmers_j, reverse_kmers_j)
        distances[i, j] = distances[j, i] = dist

    cluster_indices = []

    clustering_order = iterative_clustering(distances)

    # Print the names of the sequences in the order they were clustered
    for i in clustering_order:
        print(fasta_files[i])
else:
    # Example of loading sequences and calculating pairwise Jaccard distances
    fasta_files = ['DQ409327.1.fna', 'MW534270.1.fna', 'NC_026992.1.fna', 'NC_023541.1.fna', 'NC_024860.1.fna']
    k = 21  # Example k-mer size
    kmers_per_sequence = {}

    # Generate k-mers for each sequence
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            kmers_per_sequence[record.id] = generate_kmers(str(record.seq), k)

    # Calculate pairwise Jaccard distances
    sequence_ids = list(kmers_per_sequence.keys())
    distances = np.zeros((len(sequence_ids), len(sequence_ids)))
    for i, j in combinations(range(len(sequence_ids)), 2):
        # Unpack the tuples of forward and reverse k-mers for each sequence
        forward_kmers_i, reverse_kmers_i = kmers_per_sequence[sequence_ids[i]]
        forward_kmers_j, reverse_kmers_j = kmers_per_sequence[sequence_ids[j]]
        # Calculate the average Jaccard distance
        dist = average_jaccard(False, forward_kmers_i, reverse_kmers_i, forward_kmers_j, reverse_kmers_j)
        distances[i, j] = distances[j, i] = dist

    # distances now contains the pairwise distances between your sequences

    print ("Distances:", distances)

    # Initially, no indices are in the cluster
    cluster_indices = []

    clustering_order = iterative_clustering(distances)
    print("Clustering Order:", clustering_order)
    #print the names of the sequences in the order they were clustered
    print("Sequence names in clustering order:", [sequence_ids[i] for i in clustering_order])

