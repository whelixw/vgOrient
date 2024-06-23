import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def random_transform_fasta(input_file, output_dir, reverse, transformation_log):
    """
    Optionally reverses and rotates a fasta file.
    - Chance to reverse complement each sequence if reverse is set to True.
    - Randomly rotates the sequence by a position between 0 and the length of the sequence.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    
    transformed_records = []
    for record in records:
        new_seq = record.seq
        transformation_info = {}
        transformation_info['file'] = os.path.basename(input_file)
        transformation_info['id'] = record.id

        if reverse and random.random() < 0.5:
            new_seq = new_seq.reverse_complement()
            transformation_info['flipped'] = "Yes"
        else:
            transformation_info['flipped'] = "No"

        rotation_cut = random.randint(0, len(new_seq) - 1)
        new_seq = new_seq[rotation_cut:] + new_seq[:rotation_cut]
        transformation_info['rotation_cut'] = rotation_cut
        
        transformed_record = SeqRecord(
            new_seq,
            id=record.id,
            description=record.description
        )
        transformed_records.append(transformed_record)
        transformation_log.append(transformation_info)
    
    output_file = f"{os.path.splitext(os.path.basename(input_file))[0]}.fasta"
    output_path = os.path.join(output_dir, output_file)
    
    SeqIO.write(transformed_records, output_path, "fasta")

def process_directory(directory, output_prefix, reverse, num_replicates, mute_stdout):
    for replicate_number in range(1, num_replicates + 1):
        output_dir = os.path.join(directory, f"{output_prefix}_{replicate_number}")
        os.makedirs(output_dir, exist_ok=True)

        transformation_log = []
        for filename in os.listdir(directory):
            if filename.endswith(".fasta"):
                input_file = os.path.join(directory, filename)
                random_transform_fasta(input_file, output_dir, reverse, transformation_log)
        
        # Log all transformations to a single file per replicate
        details_file = os.path.join(output_dir, "replicate_transformations.txt")
        with open(details_file, 'w') as df:
            for detail in transformation_log:
                df.write(f"File: {detail['file']}, ID: {detail['id']}, Flipped: {detail['flipped']}, Rotation Cut: {detail['rotation_cut']}\n")
        
        if not mute_stdout:
            print(f"Completed replicate {replicate_number} in {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Processes fasta files to optionally reverse and randomly rotate sequences.")
    parser.add_argument("directory", type=str, help="Directory containing FASTA files to process.")
    parser.add_argument("--output_prefix", type=str, default="transformed", help="Prefix for the output folder.")
    parser.add_argument("--reverse", action="store_true", help="Flag to reverse complement sequences with 50%% chance.")
    parser.add_argument("-n", "--num_replicates", type=int, default=1, help="Number of replicates to generate.")
    parser.add_argument("--mute", action="store_true", help="Mute stdout for transformations details.")

    args = parser.parse_args()

    process_directory(args.directory, args.output_prefix, args.reverse, args.num_replicates, args.mute)

if __name__ == "__main__":
    main()
