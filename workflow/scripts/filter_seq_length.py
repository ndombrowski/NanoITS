import argparse
from Bio import SeqIO
import os

def filter_sequences(input_file, output_file, min_length):
    # Parse input fasta file and filter sequences
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) >= min_length:
            records.append(record)

    # Write filtered sequences to output fasta file
    with open(output_file, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    return len(records)

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Filter sequences based on minimum length.")
    parser.add_argument("-i", "--input_file", required=True, help="Input fasta file")
    parser.add_argument("-o", "--output_file", required=True, help="Output fasta file")
    parser.add_argument("-l", "--min_length", type=int, required=True, help="Minimum sequence length to keep")

    # Parse command line arguments
    args = parser.parse_args()

    # Get the number of input sequences
    print(f"Start working on: {args.input_file}")
    input_seqs_count = len(list(SeqIO.parse(args.input_file, "fasta")))
    base_file_name = os.path.basename(args.input_file)
    
    # Filter sequences and get the number of output sequences
    output_seqs_count = filter_sequences(args.input_file, args.output_file, args.min_length)

    # Print output message
    print(f"File {base_file_name} filtered successfully.\n Number input seqs: {input_seqs_count}\n Number output seqs: {output_seqs_count}.\n Output saved to {args.output_file}")

if __name__ == "__main__":
    main()
