process fastqc {
    
    input:
    tuple val(sample), path(reads)
    val(stage)

    output:
    tuple val(sample), path("stats/${stage}/fastqc/*/*fastqc_data.txt"), emit: stats
    tuple val(sample), path("stats/${stage}/read_counts/${sample.id}.${stage}.txt"), emit: counts

    script:

    def compression = ""
    def process_r2 = ""

    if (sample.is_paired) {
        compression = reads[0].endsWith(".gz") ? "gz" : "bz2"
        process_r2 = "fastqc -t $task.cpus --extract --outdir=fastqc ${sample.id}_R2.fastq.${compression} && mv fastqc/${sample.id}_R2_fastqc/fastqc_data.txt fastqc/${sample.id}_R2_fastqc/${sample.id}_R2_fastqc_data.txt"
    } else {
        compression = reads.endsWith(".gz") ? "gz" : "bz2"
    }

    """
    set -e -o pipefail
    mkdir -p stats/${stage}/read_counts fastqc/
    fastqc -t $task.cpus --extract --outdir=fastqc ${sample.id}_R1.fastq.${compression} && mv fastqc/${sample.id}_R1_fastqc/fastqc_data.txt fastqc/${sample.id}_R1_fastqc/${sample.id}_R1_fastqc_data.txt
    ${process_r2}

    grep "Total Sequences" fastqc/*/*data.txt > seqcount.txt
    echo \$(wc -l seqcount.txt)\$'\t'\$(head -n1 seqcount.txt | cut -f 2) > stats/${stage}/read_counts/${sample.id}.${stage}.txt
	mv fastqc stats/${stage}/
    """
}
