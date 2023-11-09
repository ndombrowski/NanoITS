rule run_itsx:
    input:
        "results/{project}/reads/2_qual_filtered_reads/3_fasta/{sample}.fasta"
    output:
        "results/{project}/itsx/{sample}/{sample}.full.fasta"
    params:
        threads=config["itsx"]["threads"]
    log:
        "logs/{project}_itsx_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/itsx_{sample}.txt"
    shell: """
    #generate temp folder
    mkdir -p results/{project}/itsx/{wildcards.sample}/temp
    
    echo "Running ITSx with the following command:" > {log}
    echo "ITSx -i {input} \
            --save_regions all \
            -o results/{project}/itsx/{wildcards.sample}/{wildcards.sample} \
            --cpu {params.threads} --preserve F --table T \
            --temp results/{project}/itsx/{wildcards.sample}/temp" >> {log}

    ITSx -i {input} \
            --save_regions all \
            -o results/{project}/itsx/{wildcards.sample}/{wildcards.sample} \
            --cpu {params.threads} --preserve F --table T \
            --temp results/{project}/itsx/{wildcards.sample}/temp >> {log} 2>&1

    #cleanup
    rm -r results/{project}/itsx/{wildcards.sample}/temp
    """



rule parse_itsx:
    input:
        "results/{project}/itsx/{sample}/{sample}.full.fasta"
    output:
        "results/{project}/itsx/{sample}/{sample}_{marker}_final.fasta"
    params:
        min_itsx_length = config["itsx"]["min_its_length"]
    log:
        "logs/{project}_parse_itsx_{marker}_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    #remove short sequences
    python {my_basedir}/scripts/filter_seq_length.py \
            -i "results/{project}/itsx/{wildcards.sample}/{wildcards.sample}.{wildcards.marker}.fasta" \
            -l {params.min_itsx_length} \
            -o results/{project}/itsx/{wildcards.sample}/{wildcards.sample}_{wildcards.marker}_temp.fasta \
            > {log} 2>&1
    
    #clean fasta header
    cut -f1 -d " " results/{project}/itsx/{wildcards.sample}/{wildcards.sample}_{wildcards.marker}_temp.fasta \
                > {output}
    
    #cleanup
    rm results/{project}/itsx/{wildcards.sample}/{wildcards.sample}_{wildcards.marker}_temp.fasta
    """



# rule summarize_itsx:
#     input:
#         "results/{project}/itsx/"
#     output:
#         "results/{project}/itsx/itsx_counts_summary.txt",
#         "results/{project}/itsx/itsx_seq_length.txt"
#     params:
#         min_itsx_length = config["itsx"]["min_its_length"]
#     log:
#         "logs/summarize_itsx.log"
#     conda:
#         "../envs/nanopore.yaml"
#     priority: -50
#     shell: """
#     #count nr of hits per marker
#     python {my_basedir}/scripts/itsx_get_counts.py \
#             --its_results/{project} {input} -l {params.min_itsx_length} \
#             -o {output[0]}

#     """



