import argparse
import glob
import os
import pandas as pd

def parse_nanoplot_stats(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:10]
    
    data = []
    for line in lines:
        parts = line.strip().split(':')
        if len(parts) == 2:
            value = parts[1].strip().replace(',', '')
            data.append(float(value) if '.' in value else int(value))

    return data

def merge_nanoplot_stats(input_pattern, output_file):
    file_paths = sorted(glob.glob(input_pattern))
    data_dict = {}

    for file_path in file_paths:
        folder_name = os.path.basename(os.path.dirname(file_path))
        data_dict[folder_name] = parse_nanoplot_stats(file_path)

    df = pd.DataFrame.from_dict(data_dict, orient='index', columns=[
        "Mean_read_length",
        "Mean_read_quality",
        "Median_read_length",
        "Median_read_quality",
        "Number_of_reads",
        "Read_length N50",
        "STDEV_read_length",
        "Total_bases"
    ])

    df.to_csv(output_file, sep='\t', float_format='%.1f')

def main():
    parser = argparse.ArgumentParser(description='Merge Nanoplot stats from multiple files.')
    parser.add_argument('-i', '--input', help='Input file pattern (e.g., results/quality_checks/1_raw/nanoplot/*/NanoStats.txt)', required=True)
    parser.add_argument('-o', '--output', help='Output file (e.g., output.txt)', required=True)
    args = parser.parse_args()

    merge_nanoplot_stats(args.input, args.output)

if __name__ == "__main__":
    main()
