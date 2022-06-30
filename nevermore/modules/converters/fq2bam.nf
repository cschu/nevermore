process fq2bam {
    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("out/${sample.id}.bam"), emit: reads

    script:
	def maxmem = task.memory.toGiga()
	def r2 = (sample.is_paired) ? "in2=${sample.id}_R2.fastq.gz" : ""

	"""
	mkdir -p out/
	reformat.sh -Xmx${maxmem}g in=${sample.id}_R1.fastq.gz ${r2} trimreaddescription=t out=stdout.bam | samtools addreplacerg -r "ID:A" -r "SM:${sample.id}" -o out/${sample.id}.bam -
	"""
}