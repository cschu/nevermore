process remove_host_kraken2 {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(fq)
	path(kraken_db)

    output:
    tuple val(sample), path("no_host/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

    script:
    def out_options = (sample.is_paired) ? "--paired --unclassified-out ${sample.id}#.fastq" : "--unclassified-out ${sample.id}_1.fastq"
    def move_r2 = (sample.is_paired) ? "gzip -c ${sample.id}_2.fastq > no_host/${sample.id}/${sample.id}_R2.fastq.gz" : ""

    """
    mkdir -p no_host/${sample.id}
    kraken2 --threads $task.cpus --db ${kraken_db} ${out_options} --output kraken_read_report.txt --report kraken_report.txt --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups} $fq

    gzip -c ${sample.id}_1.fastq > no_host/${sample.id}/${sample.id}_R1.fastq.gz
    ${move_r2}
    """
}
