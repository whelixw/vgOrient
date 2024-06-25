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

def save_fasta(fasta, filename):
    """Save the fetched FASTA sequence to a file."""
    with open(filename, 'w') as file:
        file.write(fasta)
    print(f"Saved FASTA for {filename}")

def file_exists(filename, directory):
    """Check if the file exists in the specified directory."""
    return os.path.isfile(os.path.join(directory, filename)) if directory else os.path.isfile(filename)

def main(args):
    # Set the directory for saving FASTA files
    directory = args.directory if args.directory else "."
    # Set email
    email = args.email if args.email else "default_email@example.com"

    for accession in args.accessions:
        filename = f"{accession}.fasta"
        full_path = os.path.join(directory, filename)

        if not file_exists(filename, directory):
            fasta = fetch_fasta_from_ncbi(accession, email)
            if not fasta:
                print(f"No data from NCBI for {accession}, trying ENA...")
                fasta = fetch_fasta_from_ena(accession)
            if fasta:
                save_fasta(fasta, full_path)
            else:
                print(f"Failed to retrieve FASTA for {accession}")
        else:
            print(f"File {filename} already exists. Skipping download.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch FASTA sequences from NCBI and ENA")
    parser.add_argument('accessions', nargs='+', help='List of accession numbers')
    parser.add_argument('--directory', '-d', help='Directory to save FASTA files', required=False)
    parser.add_argument('--email', '-e', help='Email for NCBI requests', required=False)
    args = parser.parse_args()
    main(args)
