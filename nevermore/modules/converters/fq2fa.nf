process fq2fa {
    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("out/${sample.id}*.fasta"), emit: reads

    script:
	def maxmem = task.memory.toGiga()
	def r2 = (sample.is_paired) ? "in2=${sample.id}_R2.fastq.gz out2=out/${sample.id}_R2.fasta" : ""

	"""
	mkdir -p out/
	reformat.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz out=out/${sample.id}_R1.fasta ${r2} trimreaddescription=t
	"""
}
