import argparse
import pandas as pd
from tabulate import tabulate
import subprocess

def extract_sample_id(filepath):
    # Extract sample ID from filepath
    return filepath.split('/')[-1].split('_')[0]

def process_file(file_path):
    # Read the file and extract relevant rows
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = {}
        for line in lines[1:9]:
            items = line.strip().split(':')
            if len(items) == 2:
                key = items[0].strip().replace(" ", "_")
                value = items[1].strip()
                data[key] = value
            else:
                print(f"Warning: Skipping invalid line in {file_path}: {line}")

        return data

def main():
    parser = argparse.ArgumentParser(description='Merge nanostat files by sample ID.')
    parser.add_argument('-i', '--input', nargs='+', help='Input nanostat files', required=True)
    parser.add_argument('-o', '--output_file', help='Output merged file in CSV format', required=True)
    args = parser.parse_args()

    # Process each input file
    data_dict = {}
    for input_file in args.input:
        sample_id = extract_sample_id(input_file)
        data_dict[sample_id] = process_file(input_file)

    # Create DataFrame
    df = pd.DataFrame.from_dict(data_dict, orient='index', columns=[
        "Mean_read_length",
        "Mean_read_quality",
        "Median_read_length",
        "Median_read_quality",
        "Number_of_reads",
        "Read_length_N50",
        "STDEV_read_length",
        "Total_bases"
    ])

    #df_transposed = df.T

    #write to csv
    df.to_csv(f"{args.output_file.replace('.pdf', '.txt')}", index_label = "stats")

    # Replace underscores in DataFrame index
    df.index = df.index.str.replace("_", "\\_")
    df.columns = df.columns.str.replace("_", "\\_")

    # Save DataFrame to a LaTeX table with underscores replaced
    latex_table = tabulate(df, tablefmt="latex_raw", headers="keys")
    #latex_table = tabulate(df, tablefmt="latex_raw", headers="keys", colalign=("center",) * len(df_transposed.columns), floatfmt=(".2f",) * len(df_transposed.columns), table_options=["landscape"])
    
    # Save the LaTeX table to the specified output file
    with open(args.output_file, "w") as f:
        f.write(latex_table)

    # Convert LaTeX to PDF using Pandoc
    subprocess.run(["pandoc", args.output_file, "-o", f"{args.output_file.replace('.tex', '.pdf')}"])


if __name__ == '__main__':
    main()

