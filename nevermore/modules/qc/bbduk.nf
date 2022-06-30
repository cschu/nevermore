process qc_bbduk {
	label 'bbduk'
	publishDir path: params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("${sample.id}/${sample.id}.bbduk_stats.txt")

    script:
    def maxmem = task.memory.toGiga()
    def read1 = "in1=${sample.id}_R1.fastq.gz out1=${sample.id}/${sample.id}_R1.fastq.gz"
    def read2 = sample.is_paired ? "in2=${sample.id}_R2.fastq.gz out2=${sample.id}/${sample.id}_R2.fastq.gz outs=${sample.id}/${sample.id}.orphans_R1.fastq.gz" : ""

	qc_params = params.qc_params_shotgun
	trim_params = "ref=${adapters} minlen=${params.qc_minlen}"

    """
    mkdir -p ${sample.id}
    bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${qc_params} ref=${adapters} stats=${sample.id}/${sample.id}.bbduk_stats.txt ${read1} ${read2}
    """
}
