rule run_kraken:
    input:
        query = "results/{project}/itsx/{sample}/{sample}_{marker}_final.fasta",
        silva_kraken_db = os.path.join(my_basedir, "db/silva/kraken/"),
        silva_unite_db = os.path.join(my_basedir, "db/unite/kraken/")
    output:
        out = "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.out",
        report = "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.report"
    params:
        threads = config["threads"],
        confidence_score = config["kraken2"]["confidence_score"],
        min_hit_nr = config["kraken2"]["min_hit_nr"]
    log:
        "logs/{project}_kraken2_{marker}_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/{project}_kraken2_{marker}_{sample}.txt"
    shell: """
    
    if [ "{wildcards.marker}" == "SSU" ]; then
        kraken_db="{input.silva_kraken_db}"
    elif [ "{wildcards.marker}" == "ITS1" ]; then
        kraken_db="{input.silva_unite_db}"
    else
        echo "Unsupported marker: {wildcards.marker}"
    fi

    echo "Running kraken2 with: "  > {log}
    echo "kraken2 -db $kraken_db --confidence {params.confidence_score} --minimum-hit-groups {params.min_hit_nr} --threads {params.threads}"  >> {log}
    
    echo "Starting kraken2 --version"
    kraken2_version=$(kraken2 --version)
    echo "kraken2 version is $kraken2_version" >> {log}

    kraken2 --db $kraken_db \
            --confidence {params.confidence_score} \
            --minimum-hit-groups {params.min_hit_nr} \
            --output {output.out} \
            --report {output.report} \
            --threads {params.threads}  \
            {input.query} 2>> {log}

    """


rule kraken_to_lca:
    input:
        results = "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.out",
        silva_kraken_tax = os.path.join(my_basedir, "db/silva/kraken/taxID_to_tax.txt"),
        unite_kraken_tax = os.path.join(my_basedir, "db/unite/kraken/taxID_to_tax.txt")
    output:
        "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.taxlist"
    log:
        "logs/{project}_kraken2_to_lca_{marker}_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    
    if [ "{wildcards.marker}" == "SSU" ]; then
        reference_tax="{input.silva_kraken_tax}"
    elif [ "{wildcards.marker}" == "ITS1" ]; then
        reference_tax="{input.unite_kraken_tax}"
    else
        echo "Unsupported marker: {wildcards.marker}"
        exit 1
    fi

    echo "Used database for {wildcards.marker} in kraken_to_lca: $reference_tax"  > {log}

    python3 {my_basedir}/scripts/tolca.py \
        -b <(awk -v OFS="\t" '{{print $2, $3}}' {input.results}) \
        -t  $reference_tax \
        -l {output} \
        -c 0.55 2>> {log}
    """



rule kraken_to_matrix:
    input:
        "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.taxlist"
    output:
        "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.otumat",
        "results/{project}/classification/kraken2/{sample}/{sample}_{marker}_kraken2.taxmat"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    
    python {my_basedir}/scripts/tomat.py \
            -l {input}
    """