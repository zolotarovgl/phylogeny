import sys
from Bio import SeqIO

def extract_sequences(fasta_file, ids_file, output_file):
    # Read the list of IDs from the file
    with open(ids_file, 'r') as f:
        ids_to_extract = set(line.strip() for line in f)
    
    # Parse the FASTA file and extract sequences
    with open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in ids_to_extract:
                SeqIO.write(record, output_handle, 'fasta')
    print(f'Extracted {len(ids_to_extract)} sequences from {fasta_file} to {output_file}')

if __name__ == "__main__":
    # Ensure correct number of arguments
    if len(sys.argv) != 4:
        print("Usage: python get_seq.py <fasta_file> <ids_file> <output_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    ids_file = sys.argv[2]
    output_file = sys.argv[3]
    
    extract_sequences(fasta_file, ids_file, output_file)
