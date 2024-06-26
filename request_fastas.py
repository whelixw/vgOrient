import argparse
import os
from Bio import Entrez
import requests

def fetch_fasta_from_ncbi(accession, email):
    """Fetch the FASTA from NCBI based on the accession number."""
    Entrez.email = email
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            fasta_data = handle.read()
            if fasta_data.strip() and fasta_data.startswith('>'):
                return fasta_data
    except Exception as e:
        print(f"Error fetching from NCBI for {accession}: {e}")
    return None

def fetch_fasta_from_ena(accession):
    """Fetch the FASTA from ENA based on the accession number."""
    url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
    response = requests.get(url)
    if response.ok and response.text.strip().startswith('>'):
        return response.text
    return None

def extract_name_from_header(fasta_data):
    """Extract the first word from the FASTA header (assumed to be after '>')."""
    first_line = fasta_data.partition('\n')[0]
    first_word = first_line.split()[0][1:]  # Assumes there's a word after '>'
    return first_word

def save_fasta(fasta, filename):
    """Save the fetched FASTA sequence to a file."""
    with open(filename, 'w') as file:
        file.write(fasta)
    print(f"Saved FASTA for {filename}")

def file_exists(filename, directory):
    """Check if the file exists in the specified directory."""
    return os.path.isfile(os.path.join(directory, filename)) if directory else os.path.isfile(filename)

def main(args):
    directory = args.directory if args.directory else "."
    email = args.email if args.email else "default_email@example.com"

    for accession in args.accessions:
        fasta = fetch_fasta_from_ncbi(accession, email)
        if not fasta:
            print(f"No data from NCBI for {accession}, trying ENA...")
            fasta = fetch_fasta_from_ena(accession)
        
        if fasta:
            new_filename = f"{accession}.fasta"
            if args.rename:
                new_accession = extract_name_from_header(fasta)
                new_filename = f"{new_accession}.fasta"
            full_path = os.path.join(directory, new_filename)

            if not file_exists(new_filename, directory):
                save_fasta(fasta, full_path)
            else:
                print(f"File {new_filename} already exists. Skipping download.")
        else:
            print(f"Failed to retrieve FASTA for {accession}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch FASTA sequences from NCBI and ENA")
    parser.add_argument('accessions', nargs='+', help='List of accession numbers')
    parser.add_argument('--directory', '-d', help='Directory to save FASTA files', required=False)
    parser.add_argument('--email', '-e', help='Email for NCBI requests', required=False)
    parser.add_argument('--rename', action='store_true', help='Rename output files based on the first word in the FASTA header')
    args = parser.parse_args()
    main(args)

