import subprocess
import argparse
import time
import logging
import os

def setup_logging(log_file, log_level=logging.INFO):
    """Set up the logging configuration."""
    logging.basicConfig(filename=log_file, level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser(description="Run kmer_jaccard.py and VG_diterative.py")
    parser.add_argument('input_files', nargs='+', help='Input file paths for kmer_jaccard.py')
    parser.add_argument('--output', '-o', default='output.txt', help='Output file name for kmer_jaccard.py')
    parser.add_argument('--vg_output_dir', help='Output directory for VG_diterative.py')
    parser.add_argument('--kmer_size', '-k', default='11', help='kmer size for computing jaccard similarity')
    parser.add_argument('--orientation', action='store_true', help='Use specific orientation in kmer_jaccard.py')
    parser.add_argument('--band_width', '-w', type=int, default=512, help='Band width for VG mapping.')
    parser.add_argument('--min_match_length', '-m', type=int, default=512, help='Minimum match length for VG mapping.')
    parser.add_argument('--append_wm', action='store_true', default=True, help='Append w and m values to the output directory name.')
    parser.add_argument('--log', help='Log file name for recording execution details and timings.', default='vg_wrapper.log')

    args = parser.parse_args()

    # Construct the output directory name based on user options
    base_output_dir = args.vg_output_dir if args.vg_output_dir else "vg_output"
    if args.append_wm:
        output_dir = f"{base_output_dir}_w{args.band_width}_m{args.min_match_length}"
    else:
        output_dir = base_output_dir

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set up the log file within the output directory
    log_file_path = os.path.join(output_dir, args.log)
    setup_logging(log_file_path)

    logging.info(f"Starting script with arguments: {args}")

    # Running kmer_jaccard.py
    start_time = time.time()
    kmer_cmd = ['python3', 'kmer_jaccard.py']
    if args.kmer_size:
        kmer_cmd.extend(['-k', args.kmer_size])
    if args.orientation:
        kmer_cmd.append('--orientation')
    kmer_cmd.extend(args.input_files)
    
    subprocess.run(kmer_cmd)
    kmer_duration = time.time() - start_time
    print(f"kmer_jaccard.py execution completed.")
    logging.info(f'kmer_jaccard.py completed in {kmer_duration:.2f} seconds')

    # Running VG_diterative.py
    start_time = time.time()
    vg_cmd = ['python3', 'VG_diterative.py', '-o', output_dir, args.output, '-w', str(args.band_width), '-m', str(args.min_match_length)]
    if args.append_wm:
        vg_cmd.append('--append_wm')
    subprocess.run(vg_cmd)
    vg_duration = time.time() - start_time
    print("VG_diterative.py has been executed")
    logging.info(f'VG_diterative.py completed in {vg_duration:.2f} seconds')

if __name__ == '__main__':
    main()
