rule combine_fastq:
    input:
        "input/{sample}/"
    output:
        "results/reads/1_fastq_combined/{sample}.fastq.gz",
        "logs/read_count_{sample}.txt"
    log:
        "logs/combine_fastq_{sample}.log"
    shell: """

    # Count reads in individual input files
    for file in input/{wildcards.sample}/*fastq.gz; do
        file_read_count=$(zcat $file | wc -l)
        echo "Read counts in $file: $file_read_count" >> logs/combine_fastq_{wildcards.sample}.log
    done

    # Combine fastq files
    cat input/{wildcards.sample}/*fastq.gz > {output[0]}

    # Count reads in the combined output file
    read_count=$(zcat {output[0]} | wc -l)
    echo "Read counts in {output[0]}: $read_count" >> logs/combine_fastq_{wildcards.sample}.log

    # Store read count in the variable
    echo "$read_count" > {output[1]}
    """