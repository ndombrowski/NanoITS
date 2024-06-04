def get_fastq_location(wildcards):
    return samples_table.loc[wildcards.sample, "path"]

rule run_nanostat:
    input:
        get_fastq_location
    output:
        "results/{project}/quality_checks/1_unfiltered/nanostat/{sample}_unfiltered_stats.txt"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/{project}_NanoStatUnfiltered_{sample}.txt"
    log:
        "logs/{project}_nanostat_{sample}.log"
    threads:
        config["threads"]
    shell: """
    #run tool
    NanoStat --fastq {input} --name {output} -t {threads}
    """

# rule summarize_nanostat:
#     input:
#         expand("results/{project}/quality_checks/1_unfiltered/nanostat/{sample}_unfiltered_stats.txt", sample=samples_table.index, project=config["project"])
#     output:
#         report("results/{project}/quality_checks/1_unfiltered/nanostat/merged_stats.txt",
#         caption = "../report/fig_nanostat.rst",
#         category = "Quality_unfiltered_reads"        
#         )
#     conda:
#         "../envs/nanopore.yaml"
#     log:
#         "logs/{project}_nanostat_parsing.log"
#     shell: """
#     #run tool
#     python {my_basedir}/scripts/merge_nanostats.py \
#         -i {input} -o results/{project}/quality_checks/1_unfiltered/nanostat/merged_stats.pdf
#     """


rule run_pistis:
    input:
        get_fastq_location
    output:
        report("results/{project}/quality_checks/1_unfiltered/pistis/{sample}.pdf",
        caption = "../report/fig_pistis.rst",
        category = "Quality_unfiltered_reads"
        )
    params:
        downsample = config["pistis"]["downsample"]
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["pistis"]["memory"]    
    conda:
        "../envs/nanopore.yaml"
    log:
        "logs/{project}_pistis_{sample}.log"
    benchmark:
        "results/{project}/benchmarks/{project}_PistisUnfiltered_{sample}.txt"
    shell: """
    pistis --fastq {input} --output {output} --downsample {params.downsample}
    """


# rule run_nanoplot:
#     input:
#         get_fastq_location
#     output:
#         folder = directory("results/{project}/quality_checks/1_unfiltered/nanoplot/{sample}"),
#         results = "results/{project}/quality_checks/1_unfiltered/nanoplot/{sample}/NanoStats.txt",
#         report = report(
#                 "results/{project}/quality_checks/1_unfiltered/nanoplot/{sample}/LengthvsQualityScatterPlot_dot.png", 
#                 caption="../report/fig_nanoplot_unfiltered.rst", 
#                 category="Read-processing"
#         )
#     conda:
#         "../envs/nanopore.yaml"
#     benchmark:
#         "results/{project}/benchmarks/nanoplot_unfiltered_{sample}.txt"
#     params:
#         threads = config["threads"]
#     shell: """
#     #run tool
#     NanoPlot --fastq {input} -o {output.folder} -t {params.threads} --plots dot
#     """



# rule merge_nanoplot_stats:
#     output:
#         "results/{project}/quality_checks/1_unfiltered/nanoplot/merged_stats.txt"
#     conda:
#         "../envs/nanopore.yaml"
#     priority: -50
#     shell: """
#     python3 {my_basedir}/scripts/merge_nanoplot_stats.py \
#             -i 'results/{project}/quality_checks/1_unfiltered/nanoplot/*/NanoStats.txt' \
#             -o {output}
#     """
