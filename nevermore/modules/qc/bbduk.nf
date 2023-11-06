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
    // def maxmem = task.memory.toGiga().intdiv(2) 
    def maxmem = task.memory.toGiga() 
    def compression = (reads[0].name.endsWith("gz")) ? "gz" : "bz2"

    def read2 = ""
    def orphan_check = ""

    def bb_params = params.qc.shotgun.params //.replaceAll(/maq=([0-9]+)/, "")
    
    def trim_params = "${bb_params} ref=${adapters} minlen=${params.qc.shotgun.minlen}"

    def orphan_filter = ""
    
    if (sample.is_paired) {
        def orphans = "qc_reads/${sample.id}/${sample.id}.orphans_R1.fastq.gz"
        read2 = "in2=${sample.id}_R2.fastq.${compression} out2=qc_reads/${sample.id}/${sample.id}_R2.fastq.gz outs=tmp_orphans.fq"
        orphan_filter = "bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} in=tmp_orphans.fq out=${orphans}"

        orphan_check = """
        if [[ -z "\$(gzip -dc ${orphans} | head -n 1)" ]]; then
			rm ${orphans}
		fi
        """
    }

    def read1 = "in1=${sample.id}_R1.fastq.${compression} out1=qc_reads/${sample.id}/${sample.id}_R1.fastq.gz"
    
    def stats_out = "stats=stats/qc/bbduk/${sample.id}.bbduk_stats.txt"

    """
    set -e -o pipefail

    mkdir -p qc_reads/${sample.id}/ stats/qc/bbduk/
    bbduk.sh -Xmx${maxmem}g t=${task.cpus} ${trim_params} ${stats_out} ${read1} ${read2}
    ${orphan_filter}
    ${orphan_check}

    touch qc_reads/${sample.id}/BBDUK_FINISHED
    rm -vf *.fq
    """
}
