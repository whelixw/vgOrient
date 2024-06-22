import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def random_transform_fasta(input_file, output_dir):
    """
    Randomly flips and rotates a fasta file.
    - 50% chance to reverse complement each sequence.
    - Randomly rotates the sequence by a position between 0 and the length of the sequence.
    """
    # Read the fasta file
    records = list(SeqIO.parse(input_file, "fasta"))
    
    # Apply the rotation and flip if necessary
    transformed_records = []
    for record in records:
        new_seq = record.seq
        if random.random() < 0.5:  # 50% chance to flip (reverse complement)
            new_seq = new_seq.reverse_complement()
            flip_status = "flipped"
        else:
            flip_status = "original"

        # Rotate sequence by a random cut position within its length
        rotation_cut = random.randint(0, len(new_seq) - 1)
        new_seq = new_seq[rotation_cut:] + new_seq[:rotation_cut]
        
        # Create new record with transformed sequence
        transformed_record = SeqRecord(
            new_seq,
            id=record.id,
            description=record.description
        )
        transformed_records.append(transformed_record)
    
    # Output file path with orientation and rotation in the name
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = f"{base_name}.fasta"
    output_path = os.path.join(output_dir, output_file)
    
    # Write the transformed sequences to a new file
    SeqIO.write(transformed_records, output_path, "fasta")
    print(f"Processed {input_file}: {flip_status}, rotation cut at position {rotation_cut}")

def process_directory(directory):
    # Ensure the output directory exists
    output_dir = os.path.join(directory, "transformed")
    os.makedirs(output_dir, exist_ok=True)

    # Process each fasta file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            input_file = os.path.join(directory, filename)
            random_transform_fasta(input_file, output_dir)

def main():
    parser = argparse.ArgumentParser(description="Randomly flips and rotates fasta files in a directory.")
    parser.add_argument("directory", type=str, help="Directory containing FASTA files to process.")
    
    args = parser.parse_args()
    
    process_directory(args.directory)

if __name__ == "__main__":
    main()
