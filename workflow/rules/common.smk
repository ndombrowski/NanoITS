rule download_dbs:
    output:
        silva_general_fasta = os.path.join(my_basedir,"db/silva/general/silva-ref.fasta"),
        silva_general_tax = os.path.join(my_basedir, "db/silva/general/silva-ref-taxonomy.txt"),

        unite_general_fasta = os.path.join(my_basedir,"db/unite/general/unite-ref-seqs.fna"),
        unite_general_tax = os.path.join(my_basedir, "db/unite/general/unite-ref-taxonomy.txt"),

        silva_kraken_db = directory(os.path.join(my_basedir, "db/silva/kraken/")),
        unite_kraken_db = directory(os.path.join(my_basedir, "db/unite/kraken/"))
    log:
        "logs/run_download_db.log"
    params:
        download_path = config["db_download"]
    conda:
        "../envs/nanopore.yaml"
    threads: 1
    shell: """
    wget -q -O db.tar.gz {params.download_path}
    tar -xzf db.tar.gz -C {my_basedir}/
    rm db.tar.gz 
    """


rule summarize_itsx_length:
    input:
        expand("results/{project}/itsx/{sample}/{sample}_{marker}_final.fasta",
        project=config["project"], marker=config["markers"], sample=samples_table.index)
    output:
        report("results/{project}/itsx/avg_seq_length.txt",
        caption="../report/fig_itsx.rst",
        category="ITSx_statistics"
        )
    conda:
        "../envs/nanopore.yaml"
    shell: """
    python3  {my_basedir}/scripts/calc_itsx_seq_length.py \
        --its_results results/{project}/itsx  \
        -o results/run_v1/itsx/avg_seq_length.txt
    """


rule summarize_itsx_nr:
    input:
        expand("results/{project}/itsx/{sample}/{sample}_{marker}_final.fasta",
        project=config["project"], marker=config["markers"], sample=samples_table.index)
    output:
        report("results/{project}/itsx/avg_seq_nr.txt",
        caption="../report/fig_itsx.rst",
        category="ITSx_statistics"
        )
    conda:
        "../envs/nanopore.yaml"
    shell: """
    python3  {my_basedir}/scripts/itsx_get_counts.py \
        --its_results results/{project}/itsx  \
        -o {output}
    """


rule merge_benchmarks:
    input:
        expand("results/{project}/plotting/{marker}/{marker}_Barplot_{type}_{rank}_rank.pdf", rank=config["tax_ranks"], 
            project=config["project"], marker=config["markers"], type=["counts", "ra"])
    output:
        report("results/{project}/benchmarks/merged.pdf",
        caption="../report/fig_benchmarks.rst",
        category="Run_statistcs"
        )
    conda:
        "../envs/nanopore.yaml"
    shell: """
    python {my_basedir}/scripts/merge_benchmarks.py \
        -i "results/run_v1/benchmarks/*txt" \
        -o results/run_v1/benchmarks/merged_output.txt \
        --barplots {output}
    """