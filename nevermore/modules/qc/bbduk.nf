process qc_bbduk {
	label 'bbduk'
	// publishDir path: params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("stats/qc/bbduk/${sample.id}.bbduk_stats.txt")

    script:
    def maxmem = task.memory.toGiga()
    

    def compression = ""
    def read2 = ""

    if (sample.is_paired) {
        compression = reads[0].endsWith(".gz") ? "gz" : "bz2"
        read2 = sample.is_paired ? "in2=${sample.id}_R2.fastq.gz out2=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz outs=qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz" : ""
    } else {
        compression = reads.endsWith(".gz") ? "gz" : "bz2"
    }

    def read1 = "in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
    
    trim_params = params.qc_params_shotgun + " ref=${adapters} minlen=${params.qc_minlen}"

    """
    mkdir -p qc_reads/${sample.id}
	mkdir -p stats/qc/bbduk/
    bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} stats=stats/qc/bbduk/${sample.id}.bbduk_stats.txt ${read1} ${read2}
    """
}
