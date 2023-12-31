
## NanoITS

NanoITS is a classifier for long-read Oxford Nanopore data of the eukaryotic 18S/SSU-ITS1 operon. 

When giving the tool some nanopore long-read data it will:

- Provide a quality report of the raw reads
- Check the reads for adaptors and barcodes and if present trim the reads
- Remove low-quality and short reads
- Provide a quality report for the cleaned reads
- Identify and separate both the ITS1 and 18S rRNA gene
- Classify the SSU and/or ITS1 gene using kraken2 and/or minimap2
- Generate taxonomic barplots and OTU tables

We plan to in the future extend this workflow to also include the 28S/LSU.

For a more detailed explanation, check out [the manual](https://ndombrowski.github.io/NanoITS/). 

Below you can find the full workflow:

![](img/visualization.png)


## Quick start

To run NanoITs, install conda and snakemake and clone the directory from github via:

```{python}
git clone https://github.com/ndombrowski/NanoITS.git
```

Provide your sample names and path to the samples as a comma-separated file, for example, a file looking similar as the one provided in `example_files/mapping.csv`. Sample names should be unique and consist of letter, numbers and `-` only. The barcode column can be left empty as it is not yet implemented. The path should contain the path to your demultiplexed, compressed fastq file.

Adjust `config/config.yaml` to configure the location of your mapping file as well as specify the parameters used by NanoITs.

NanoITs can then be run with:


```{python}
snakemake --use-conda --cores <nr_cores> \
  -s <path_to_NanoITS_install>/workflow/Snakefile \
  --configfile config/config.yaml \
  --conda-prefix <path_to_NanoITS_install>/workflow/.snakemake/conda  \
  --rerun-incomplete --nolock -np 
```

For a more detailed explanation, check out [the manual](https://ndombrowski.github.io/NanoITS/).