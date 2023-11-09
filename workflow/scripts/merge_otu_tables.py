import os
import pandas as pd
import argparse

def extract_barcodes_from_folder(results_folder):
    barcodes = []

    # Walk through the directory structure
    for root, dirs, files in os.walk(results_folder):
        # Extract variables from the folder paths
        for dir in dirs:
            barcode = os.path.join(root, dir).split(os.path.sep)[4]
            barcodes.append(barcode)

    return barcodes

def merge_files(barcodes, input_dir, output_file, marker, method):

    # Initialize an empty dataframe to store the merged data
    merged_df = pd.DataFrame(columns=["#NAME"])

    # Loop through each barcode
    for barcode in barcodes:
        # Define the path to the result file for the current barcode
        result_file = os.path.join(input_dir, f"{barcode}/{barcode}_{marker}_{method}.otumat")

        # Read the data from the result file into a dataframe
        current_df = pd.read_csv(result_file, sep='\t')
        
        #change header
        current_df.rename(columns={'#taxid': '#NAME', f'{method}_{barcode}': f'{barcode}'}, inplace=True)

        # Merge the data into the main dataframe based on the "#NAME" column
        merged_df = pd.merge(merged_df, current_df, on="#NAME", how="outer")
        
    # Fill NaN values with 0
    merged_df = merged_df.fillna(0)

    # Save the merged dataframe to the output file
    merged_df.to_csv(output_file, index=False, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge files with taxids and counts.")
    parser.add_argument("-i", "--input_dir", help="Input directory containing result files")
    parser.add_argument("-m", "--marker", help="Marker to analyse")
    parser.add_argument("-o", "--output_file", help="Output file to store the merged data")
    parser.add_argument('--method', required=True, help='Classification method')
    args = parser.parse_args()

    # Extract barcodes from folder structure
    barcodes = extract_barcodes_from_folder(args.input_dir)

    merge_files(barcodes, args.input_dir, args.output_file, args.marker, args.method)


