process qc_bbduk {
	label 'bbduk'

    input:
    tuple val(sample), path(reads)
	path(adapters)

    output:
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
    tuple val(sample), path("qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"), emit: orphans, optional: true
    path("stats/qc/bbduk/${sample.id}.bbduk_stats.txt")
    tuple val(sample), path("qc_reads/${sample.id}/BBDUK_FINISHED"), emit: sentinel

    script:
    def maxmem = task.memory.toGiga()
    def compression = (reads[0].name.endsWith("gz")) ? "gz" : "bz2"

    def read2 = ""
    def orphan_check = ""
    
    if (sample.is_paired) {
        def orphans = "qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"
        read2 = "in2=${sample.id}_R2.fastq.${compression} out2=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz outs=${orphans}"
        orphan_check = """
        if [[ -z "\$(gzip -dc ${orphans} | head -n 1)" ]]; then
			rm ${orphans}
		fi
        """
    }

    def read1 = "in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
    
    def trim_params = params.qc_params_shotgun + " ref=${adapters} minlen=${params.qc_minlen}"
    def stats_out = "stats=stats/qc/bbduk/${sample.id}.bbduk_stats.txt"

    """
    set -e -o pipefail

    mkdir -p qc_reads/${sample.id}/ stats/qc/bbduk/
    bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} ${stats_out} ${read1} ${read2}
    ${orphan_check}

    touch qc_reads/${sample.id}/BBDUK_FINISHED
    """
}
