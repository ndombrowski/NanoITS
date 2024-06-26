rule run_minimap:
    input:
        query = "results/{project}/itsx/{sample}/{sample}_{marker}_final.fasta",
        silva_general_fasta = os.path.join(my_basedir, "db/silva/general/silva-ref.fasta"),
        unite_general_fasta = os.path.join(my_basedir, "db/unite/general/unite-ref-seqs.fna")
    output:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.paf"
    threads:
        config["minimap2"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["minimap2"]["memory"]
    log:
        "logs/{project}_minimap2_{marker}_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/minimap2_{marker}_{sample}.txt"
    shell: """

    if [ "{wildcards.marker}" == "SSU" ]; then
        reference_fasta="{input.silva_general_fasta}"
    elif [ "{wildcards.marker}" == "ITS1" ] || [ "{wildcards.marker}" == "ITS2" ]; then
        reference_fasta="{input.unite_general_fasta}"
    else
        echo "Unsupported marker: {wildcards.marker}"
        exit
    fi

    echo "Running minimap2 with:"  >> {log}
    echo "minimap2 -cx map-ont -t {threads} -N 10 -K 25M $reference_fasta" >> {log}
  
    echo "Starting minimap2"  >> {log}
    minimap_version=$(minimap2 --version)
    echo "minimap version is $minimap_version" >> {log}

    minimap2 -cx map-ont -t {threads} \
            -N 10 -K 25M \
            $reference_fasta \
            {input.query} \
            -o {output} &>> {log}

    """



rule parse_minimap:
    input:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.paf"
    output:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.out"
    log:
        "logs/{project}_parse_minimap_{marker}_{sample}.log"
    params:
        coverage_threshold=config["minimap2"]["coverage_threshold"]
    conda:
        "../envs/nanopore.yaml"
    shell: """
    sed 's/AS:i://' {input} | \
        awk -F'\t' -v OFS='\t' '{{len=$2; start=$3; end=$4; print $1, $6, $12, $15, len, start, end, $10, $11, (end-start+1)/len*100, $10/$11}}' | \
        python3 {my_basedir}/scripts/filter_paf.py -i - -o - -c {params.coverage_threshold} | \
        awk -v OFS="\t" '{{print $1, $2}}' >  {output}
    """



rule minimap_to_lca:
    input:
        results = "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.out",
        silva_general_tax = os.path.join(my_basedir, "db/silva/general/silva-ref-taxonomy.txt"),
        unite_general_tax = os.path.join(my_basedir, "db/unite/general/unite-ref-taxonomy.txt")
    output:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.taxlist"
    log:
        "logs/{project}_minimap2_to_lca_{marker}_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    
    if [ "{wildcards.marker}" == "SSU" ]; then
        reference_tax="{input.silva_general_tax}"
    elif [ "{wildcards.marker}" == "ITS1" ] || [ "{wildcards.marker}" == "ITS2" ]; then
        reference_tax="{input.unite_general_tax}"
    else
        echo "Unsupported marker: {wildcards.marker}"
    fi

    echo "Used database for {wildcards.marker} in minimap2 to lca: $reference_tax" > {log}

    python {my_basedir}/scripts/tolca.py \
            -b {input.results} \
            -t $reference_tax \
            -l {output} \
            -c 0.55
    """



rule minimap_to_matrix:
    input:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.taxlist"
    output:
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.otumat",
        "results/{project}/classification/minimap2/{sample}/{sample}_{marker}_minimap2.taxmat"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    
    python {my_basedir}/scripts/tomat.py \
            -l {input}
    """