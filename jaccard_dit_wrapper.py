import subprocess
import argparse
import time
import logging

def main():
    parser = argparse.ArgumentParser(description="Run kmer_jaccard.py and VG_diterative.py")
    parser.add_argument('input_files', nargs='+', help='Input file paths for kmer_jaccard.py')
    parser.add_argument('--output', '-o', default='output.txt', help='Output file name for kmer_jaccard.py')
    parser.add_argument('--vg_output_dir', help='Output directory for VG_diterative.py')
    parser.add_argument('--kmer_size', '-k', default='11', help='kmer size for computing jaccard similarity')
    parser.add_argument('--log', help='Log file to record execution details and timings')

    args = parser.parse_args()

    # Set up logging
    if args.log:
        logging.basicConfig(filename=args.log, level=logging.INFO)
        logging.info(f'Starting script with arguments: {args}')

    # Running kmer_jaccard.py
    start_time = time.time()
    if args.kmer_size:
        kmer_cmd = ['python3', 'kmer_jaccard.py', '-k', args.kmer_size] + args.input_files
    else:
        kmer_cmd = ['python3', 'kmer_jaccard.py'] + args.input_files
    with open(args.output, 'w') as f:
        subprocess.run(kmer_cmd, stdout=f)
    kmer_duration = time.time() - start_time
    print(f"kmer_jaccard.py output saved to {args.output}")
    if args.log:
        logging.info(f'kmer_jaccard.py completed in {kmer_duration:.2f} seconds')

    # Running VG_diterative.py
    start_time = time.time()
    vg_cmd = ['python3', 'VG_diterative.py', '-o', args.vg_output_dir, args.output] if args.vg_output_dir else ['python3', 'VG_diterative.py', args.output]
    subprocess.run(vg_cmd)
    vg_duration = time.time() - start_time
    print("VG_diterative.py has been executed")
    if args.log:
        logging.info(f'VG_diterative.py completed in {vg_duration:.2f} seconds')

if __name__ == '__main__':
    main()
