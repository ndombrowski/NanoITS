import os
import argparse

def extract_barcodes_from_folder(itsx_results_folder):
    barcodes = []

    # Walk through the directory structure
    for root, dirs, files in os.walk(itsx_results_folder):
        # Extract barcodes from the third folder in the path
        for dir in dirs:
            barcode = os.path.join(root, dir).split(os.path.sep)[3]
            barcodes.append(barcode)

    return barcodes


def get_read_count(filename):
    count = 0
    if os.path.exists(filename):
        # Use grep to count occurrences of ">"
        command = f'grep -c ">" {filename}'
        result = os.popen(command).read()
        count = int(result.strip())
    return count


def generate_table(barcodes, itsx_results_folder, out_file):
    # Initialize the table headers
    table = ['BC\tITS1_before\tITS1_after\tITS2_before\tITS2_after\tLSU_before\tLSU_after\tSSU_before\tSSU_after\t5_8S_before\t5_8S_after']

    # Loop through each barcode
    for barcode in barcodes:
        # Initialize the row with the barcode
        row = [barcode]

        # Loop through each sequence type
        for seq_type in ['ITS1', 'ITS2', 'LSU', 'SSU', '5_8S']:
            # Construct the filename
            before_filename = os.path.join(itsx_results_folder, barcode, f'{barcode}.{seq_type}.fasta')
            after_filename = os.path.join(itsx_results_folder, barcode, f'{barcode}_{seq_type}_final.fasta')

            # Read the count before filtering using grep
            count_before = get_read_count(before_filename)

            # Read the count after filtering using grep
            count_after = get_read_count(after_filename)

            # Add the counts to the row
            row.extend([str(count_before), str(count_after)])

        # Add the row to the table
        table.append('\t'.join(row))

    # Sort the table by barcode
    sorted_table = sorted(table)

    # Write the table to the output file
    with open(out_file, 'w') as output:
        output.write('\n'.join(sorted_table))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a table of sequence counts for different barcodes. Notice: Expect files to be in format barcode_markerID_final.fasta')
    parser.add_argument('--its_results', required=True, help='Folder containing ITS results subfolders')
    parser.add_argument('-o', '--out_file', required=True, help='Output file for the table')

    args = parser.parse_args()

    # Extract barcodes from folder structure
    barcodes = extract_barcodes_from_folder(args.its_results)

    # Generate table using extracted barcodes
    generate_table(barcodes, args.its_results, args.out_file)
