import pandas as pd
import glob
import argparse
import matplotlib
matplotlib.use('Agg')  # Use the 'Agg' backend which does not require a display
import matplotlib.pyplot as plt
import seaborn as sns

def extract_sample_id(filepath):
    # Extract sample ID from filepath
    return filepath.split('/')[-1].split('.txt')[0]


def process_file(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True)
    return df


def create_plots2(subset_df, output_path=None):
    # Set the style for seaborn
    #set sns style
    sns.set(style="whitegrid")

    #manage the plot layout
    g = sns.FacetGrid(subset_df, row="variable", aspect = 3, sharey=False)

    #add the data
    g.map_dataframe(sns.stripplot, "process", "value", alpha = 0.3, color="#338844",  edgecolor="white")

    # Set titles and labels
    g.set_xticklabels(rotation=45)
    g.set_axis_labels('Job ID', '')

    #change labels
    axes = g.axes.flatten()
    axes[0].set_ylabel("Max RSS")
    axes[1].set_ylabel("Mean Load")
    axes[2].set_ylabel("Running time (min)")
    g.set_titles(row_template="")

    g.fig.tight_layout()
    
    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(output_path)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Merge benchmark files.')
    parser.add_argument('-i', '--input', help='Input file pattern (e.g., "results/run_v1/benchmarks/*txt")', required=True)
    parser.add_argument('-o', '--output', help='Output file path (e.g., "merged_output.txt")', required=True)
    parser.add_argument('--barplots', help='Output path for barplots', default=None)
    args = parser.parse_args()

    # Get file paths based on the input pattern
    filepath = glob.glob(args.input)

    data_dict = {}

    for input_file in filepath:
        sample_id = extract_sample_id(input_file)
        data_dict[sample_id] = process_file(input_file)

    # Combine the dataframes
    merged_df = pd.concat(data_dict.values(), keys=data_dict.keys())

    # Reset the index and rename columns
    merged_df.reset_index(level=0, inplace=True)
    merged_df.reset_index(drop=True, inplace=True)
    merged_df.rename(columns={'level_0': 'job_id'}, inplace=True)

    #calculate runtime in min
    merged_df["min"]=merged_df["s"]/60
    merged_df['sample_id'] = merged_df['job_id'].str.split('_').str[-1]
    merged_df['process'] = merged_df['job_id'].str.rsplit('_', 1).str[0]

    # Sort DataFrame by job_id
    merged_df = merged_df.sort_values('job_id')
    
    # Reshape the DataFrame to long-form
    merged_df_long = pd.melt(merged_df, id_vars=['job_id', 'sample_id', 'process'], var_name='variable', value_name='value')
    merged_df_long['process'] = merged_df_long['process'].astype('category')

    #subset for vars of interest
    subset_df = merged_df_long[merged_df_long['variable'].isin(['min', 'mean_load', 'max_rss'])]

    # Save to the output file
    merged_df.to_csv(args.output, index=False, sep=',')

    # Create and display or save plots
    create_plots2(subset_df, args.barplots)

if __name__ == "__main__":
    main()


