import os
import argparse
from Bio import SeqIO


def extract_barcodes_from_folder(itsx_results_folder):
    barcodes = []

    # Walk through the directory structure
    for root, dirs, files in os.walk(itsx_results_folder):
        # Extract barcodes from the third folder in the path
        for dir in dirs:
            barcode = os.path.join(root, dir).split(os.path.sep)[3]
            barcodes.append(barcode)

    return barcodes


def calculate_average_length(barcodes, itsx_results_folder, out_file):

    # Initialize the table headers
    table = ['BC\tITS1_before\tITS1_after\tITS2_before\tITS2_after\tLSU_before\tLSU_after\tSSU_before\tSSU_after\t5_8S_before\t5_8S_after']

    # Loop through each barcode
    for barcode in barcodes:
        # Initialize the row with the barcode
        row = [barcode]

        # Loop through each sequence type
        for seq_type in ['ITS1', 'ITS2', 'LSU', 'SSU', '5_8S']:
            # Construct the filenames for before and after filtering
            before_filename = os.path.join(itsx_results_folder, barcode, f'{barcode}.{seq_type}.fasta')
            after_filename = os.path.join(itsx_results_folder, barcode, f'{barcode}_{seq_type}_final.fasta')

            # Read and calculate the average sequence length before filtering using Biopython
            avg_length_before = get_average_length(before_filename)

            # Read and calculate the average sequence length after filtering using Biopython
            avg_length_after = get_average_length(after_filename)

            # Add the average lengths to the row
            row.extend([str(round(avg_length_before, 2)), str(round(avg_length_after, 2))])

        # Add the row to the table
        table.append('\t'.join(row))

    # Sort the table by barcode
    sorted_table = sorted(table)

    # Write the table to the output file
    with open(out_file, 'w') as output:
        output.write('\n'.join(sorted_table))

def get_average_length(filename):
 
    avg_length = 0
    count = 0
    if os.path.exists(filename):
        records = SeqIO.parse(filename, 'fasta')
        for record in records:
            avg_length += len(record.seq)
            count += 1

        if count > 0:
            avg_length /= count

    return avg_length

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a table of average sequence lengths for different barcodes.')
    parser.add_argument('--its_results', required=True, help='Folder containing ITS results subfolders')
    parser.add_argument('-o', '--out_file', required=True, help='Output file for the table')

    args = parser.parse_args()
    
    barcodes = extract_barcodes_from_folder(args.its_results)
    calculate_average_length(barcodes, args.its_results, args.out_file)
