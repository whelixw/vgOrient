from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import random
import argparse
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description="Make replicate sequences", epilog=":)")

parser.add_argument("--input", type=str, help="Name of FASTA input file", required=True)
parser.add_argument("--output", type=str, help="Name of FASTA output file", required=True)
parser.add_argument("--log", type=str, help="Name of logfile", required=True)
parser.add_argument("--revcomp", type=bool, help="To reverse complement or not (50% chance)")
parser.add_argument("--shift", type=bool, help="To shift or not")

args = parser.parse_args()

input = args.input
output = args.output
logfile = args.log
revcomp = args.revcomp
shift = args.shift


write_handle = open(output, "w")
writer = FastaWriter(write_handle, wrap=60)

with open(logfile, "w") as log:

    with open(input, "r") as handle:

        for index, record in enumerate(SeqIO.parse(handle, "fasta")):

            # Reverse complement
            if revcomp:
                if random.random() < 0.5:
                    record.seq = record.seq.reverse_complement()
                    print(record.id, "\t", "Reverse complemented", file=log)

            # Introduce random shift in sequence
            if shift:
                shift = random.randint(0, len(record.seq))
                left = record.seq[:shift]
                right = record.seq[shift:]
                record.seq = right + left
                print(record.id, "\t", "Shifted", "\t", shift, file=log)

            writer.write_record(record)


    write_handle.close()