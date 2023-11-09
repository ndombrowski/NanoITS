rule merge_tables:
    input:
        expand("results/{project}/classification/{classifier}/{sample}/{sample}_{marker}_{classifier}.otumat",
            sample=samples_table.index, marker=config["markers"], classifier=config["classifiers"], project=config["project"])
    output:
        "results/{project}/classification/{classifier}/{marker}.merged.outmat.tsv"
    conda:
        "../envs/nanopore.yaml"
    threads: 1
    shell: """
    python3 {my_basedir}/scripts/merge_otu_tables.py \
            -i results/{wildcards.project}/classification/{wildcards.classifier} \
            -o {output} \
            -m {wildcards.marker} \
            --method {wildcards.classifier}
    """

markers = config["markers"]

rule generate_barplots:
    input:
        expand("results/{project}/classification/{classifier}/{marker}.merged.outmat.tsv", 
                classifier=config["classifiers"], project=config["project"], marker=config["markers"])
    output:
        report(expand("results/{project}/plotting/{marker}/{marker}_Barplot_{type}_{rank}_rank.pdf", 
                rank=config["tax_ranks"], project=config["project"], marker=config["markers"], type=["counts", "ra"]),
                caption="../report/fig_barplot.rst",
                category="Classification"
        )
    conda:
        "../envs/r_for_amplicon.yaml"
    threads: 1
    shell: """
    for i in {markers}; do
        echo $i
        mkdir -p results/{project}/tables
        Rscript {my_basedir}/scripts/otu_analysis.R results/{project}/classification/*/${{i}}.merged.outmat.tsv ${{i}} results/{project}
    done
    """