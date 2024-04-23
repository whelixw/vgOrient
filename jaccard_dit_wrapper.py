import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run kmer_jaccard.py and VG_diterative.py")
    parser.add_argument('input_files', nargs='+', help='Input file paths for kmer_jaccard.py')
    parser.add_argument('--output', '-o', default='output.txt', help='Output file name for kmer_jaccard.py')
    parser.add_argument('--vg_output_dir', help='Output directory for VG_diterative.py')
    parser.add_argument('--kmer_size', '-k', default='11', help='kmer size for computing jaccard similarity')
    
    args = parser.parse_args()

    # Running kmer_jaccard.py
    #kmer_cmd = ['python3', 'kmer_jaccard.py'] + args.input_files
    if args.kmer_size:
        kmer_cmd = ['python3', 'kmer_jaccard.py', '-k', args.kmer_size] + args.input_files
    else:
        kmer_cmd = ['python3', 'kmer_jaccard.py'] + args.input_files
    with open(args.output, 'w') as f:
        subprocess.run(kmer_cmd, stdout=f)
    print(f"kmer_jaccard.py output saved to {args.output}")

    # Running VG_diterative.py
    vg_cmd = ['python3', 'VG_diterative.py', '-o', args.vg_output_dir, args.output] if args.vg_output_dir else ['python3', 'VG_diterative.py', args.output]
    subprocess.run(vg_cmd)
    print("VG_diterative.py has been executed")

if __name__ == '__main__':
    main()
