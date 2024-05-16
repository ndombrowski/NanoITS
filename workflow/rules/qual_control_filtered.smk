rule run_nanostat_filtered:
    input:
        "results/{project}/reads/2_qual_filtered_reads/2_chopper/{sample}.fastq.gz"
    output:
        "results/{project}/quality_checks/2_filtered/nanostat/{sample}_unfiltered_stats.txt"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/{project}_NanoStatFiltered_{sample}.txt"
    params:
        threads = config["threads"]
    shell: """
    #run tool
    NanoStat --fastq {input} --name {output} -t {params.threads}
    """


# rule summarize_nanostat_filtered:
#     input:
#         expand("results/{project}/quality_checks/2_filtered/nanostat/{sample}_unfiltered_stats.txt", sample=samples_table.index, project=config["project"])
#     output:
#         report("results/{project}/quality_checks/2_filtered/nanostat/merged_stats.txt",
#         caption = "../report/fig_nanostat.rst",
#         category = "Quality_filtered_reads"        
#         )
#     conda:
#         "../envs/nanopore.yaml"
#     shell: """
#     #run tool
#     python {my_basedir}/scripts/merge_nanostats.py \
#         -i {input} -o results/{project}/quality_checks/2_filtered/nanostat/merged_stats.pdf
#     """


rule run_pistis_filtered:
    input:
        "results/{project}/reads/2_qual_filtered_reads/2_chopper/{sample}.fastq.gz"
    output:
        report("results/{project}/quality_checks/2_filtered/pistis/{sample}.pdf",
        caption = "../report/fig_pistis.rst",
        category = "Quality_filtered_reads"
        )
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/{project}_PistisFiltered_{sample}.txt"
    params:
        downsample = config["pistis"]["downsample"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["pistis"]["memory"]    

    shell: """
    pistis --fastq {input} --output {output} --downsample {params.downsample}
    """


# rule run_nanoplot_filtered:
#     input:
#         "results/{project}/reads/2_qual_filtered_reads/2_chopper/{sample}.fastq.gz"
#     output:
#         folder = directory("results/{project}/quality_checks/2_filtered/nanoplot/{sample}"),
#         results = "results/{project}/quality_checks/2_filtered/nanoplot/{sample}/NanoStats.txt",
#         report = report(
#                 "results/{project}/quality_checks/2_filtered/nanoplot/{sample}/LengthvsQualityScatterPlot_dot.png", 
#                 caption="../report/fig_nanoplot_unfiltered.rst", 
#                 category="Read-processing"
#         )
#     conda:
#         "../envs/nanopore.yaml"
#     benchmark:
#         "results/{project}/benchmarks/nanoplot_filtered_{sample}.txt"
#     params:
#         threads = config["threads"]
#     shell: """
#     #run tool
#     NanoPlot --fastq {input} -o {output.folder} -t {params.threads} --plots dot
#     """


# rule merge_nanoplot_stats_filtered:
#     output:
#         "results/{project}/quality_checks/2_filtered/nanoplot/merged_stats.txt"
#     conda:
#         "../envs/nanopore.yaml"
#     priority: -50
#     shell: """
#     python3 {my_basedir}/scripts/merge_nanoplot_stats.py \
#             -i 'results/{project}/quality_checks/2_filtered/nanoplot/*/NanoStats.txt' \
#             -o {output}
#     """