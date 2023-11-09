def get_fastq_location(wildcards):
    return samples_table.loc[wildcards.sample, "path"]


rule run_porechop:
    input:
        get_fastq_location
    output:
        "results/{project}/reads/2_qual_filtered_reads/1_porechop/{sample}.fastq.gz"
    log:
        "logs/{project}_porechop_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    benchmark:
        "results/{project}/benchmarks/{project}_porechop_{sample}.txt"
    params:
        threads = config["filtering"]["threads"]
    shell: """
    echo "Running porechop with the following command:" > {log}
    echo "porechop --input {input} --output {output} --threads {params.threads} --discard_middle" >> {log}

    porechop --input {input} \
            --output {output} \
            --threads {params.threads} \
            --discard_middle \
            > {log} 2>&1
    """


rule run_chopper:
    input:
        "results/{project}/reads/2_qual_filtered_reads/1_porechop/{sample}.fastq.gz"
    output:
        "results/{project}/reads/2_qual_filtered_reads/2_chopper/{sample}.fastq.gz"
    params:
        threads = config["filtering"]["threads"],
        quality_score = config["filtering"]["quality_score"],
        min_length = config["filtering"]["min_length"],
        max_length = config["filtering"]["max_length"],
        headcrop = config["filtering"]["headcrop"],
        tailcrop = config["filtering"]["tailcrop"]
    log:
        "logs/{project}_chopper_{sample}.log"
    conda:
        "../envs/nanopore.yaml"
    shell: """
    echo "Running chopper with the following command:" > {log}
    echo "gunzip -c {input} | chopper -q {params.quality_score} -l {params.min_length} --maxlength {params.max_length}  --threads {params.threads}" >> {log}

    gunzip -c {input} |\
            chopper -q {params.quality_score} \
            --headcrop {params.headcrop} --tailcrop {params.tailcrop}   \
            -l {params.min_length} --maxlength {params.max_length}  --threads {params.threads} |\
            gzip > {output} 
    """


rule generate_fasta:
    input:
        "results/{project}/reads/2_qual_filtered_reads/2_chopper/{sample}.fastq.gz"
    output:
        "results/{project}/reads/2_qual_filtered_reads/3_fasta/{sample}.fasta"
    shell: """
    zcat < {input} | \
            sed -n '1~4s/^@/>/p;2~4p' \
            > {output}
    """
