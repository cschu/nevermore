#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (!params.publish_mode) {
	params.publish_mode = "link"
}

if (!params.output_dir) {
	params.output_dir = "nevermore_out"
}

output_dir = "${params.output_dir}"


process prepare_fastqs {
	input:
	tuple val(sample), path(fq)

	output:
	tuple val(sample), path("out/${sample}_R*.fastq.gz"), emit: reads

	script:
	if (fq.size() == 2) {
		"""
		mkdir -p out
		ln -sf ../${fq[0]} out/${sample}_R1.fastq.gz
		ln -sf ../${fq[1]} out/${sample}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p out
		ln -sf ../${fq[0]} out/${sample}_R1.fastq.gz
		"""
	}
}


process qc_preprocess {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val("${sample}"), path("${sample}/${sample}.qc_*.fastq.gz"), emit: qc_reads

	script:
	def qc_params = "qtrim=rl trimq=25 maq=25 minlen=45"


	if (reads.size() == 2) {
		"""
		maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
		mkdir -p ${sample}

		bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} in=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz out=${sample}/${sample}.qc_R1.fastq.gz out2=${sample}/${sample}.qc_R2.fastq.gz outs=${sample}/${sample}.qc_O.fastq.gz
		"""
	} else {
		"""
		maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
		mkdir -p ${sample}

		bbduk.sh -Xmx\$maxmem t=$task.cpus ${qc_params} in=${sample}_R1.fastq.gz out=${sample}/${sample}.qc_U.fastq.gz
		"""
	}
}


process qc_bbmerge {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}*.fastq.gz"), emit: merged_reads

	script:
	def merge_params = "rsem=t extend2=20 iterations=5 ecct vstrict"
	"""
	maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
	mkdir -p ${sample}

	bbmerge.sh -Xmx\$maxmem t=$task.cpus ${merge_params} in=${sample}.qc_R1.fastq.gz in2=${sample}.qc_R2.fastq.gz out=${sample}/${sample}.merged_M.fastq.gz outu1=${sample}/${sample}.merged_R1.fastq.gz outu2=${sample}/${sample}.merged_R2.fastq.gz
	"""
}


process concat_singles {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}_concat_singles.fastq.gz"), emit: concat_reads

	script:
	"""
	mkdir -p $sample
	cat ${reads} > ${sample}/${sample}_concat_singles.fastq.gz
	"""
}


process decontaminate {
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample}/${sample}*.fastq.gz"), emit: decon_reads

	script:
	if (reads.size() == 2) {
		"""
		cpus=\$(expr \"$task.cpus\" - 4)
		mkdir -p $sample

		bwa mem -t \$cpus ${params.human_ref} ${sample}.merged_R1.fastq.gz ${sample}.merged_R2.fastq.gz | samtools collate -@ 2 -f -O - | samtools fastq -f 4 -0 ${sample}_other.fastq.gz -1 ${sample}.decon_R1.fastq.gz -2 ${sample}.decon_R2.fastq.gz
		"""
	} else {
		"""
		cpus=\$(expr \"$task.cpus\" - 4)
		mkdir -p $sample

		bwa mem -t \$cpus ${params.human_ref} ${sample}_concat_singles.fastq.gz | samtools collate -@ 2 -f -O - | samtools fastq -f 4 -s ${sample}_other.fastq.gz
		"""
	}
}

/*

process decontaminate {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads1)
	path(reads2)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.no_human_R1.fastq.gz", emit: r1_fq
	path "${sample}/${sample}.no_human_R2.fastq.gz", optional: true, emit: r2_fq

	script:
	def r1 = reads1
	def r2 = file(reads2) ? reads2 : ""

	def r1_out = "-1 ${sample}/${sample}.no_human_R1.fastq.gz"
	def r2_out = file(reads2) ? "-2 ${sample}/${sample}.no_human_R2.fastq.gz" : ""

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -t \$cpus ${params.human_ref} ${r1} ${r2} | samtools collate --threads 2 -f -O - | samtools fastq -@ 2 -f 4 ${r1_out} ${r2_out} -s singletons.fastq.gz -
	"""
}

process decontaminate_singles {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.no_human_U.fastq.gz", emit: u_fq

	script:

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	mkdir -p $sample
	bwa mem -t \$cpus ${params.human_ref} ${reads} | samtools collate --threads 2 -f -O - | samtools fastq -@ 2 -f 4 -s ${sample}/${sample}.no_human_U.fastq.gz -
	"""

}

process align {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads1)
	path(reads2)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}.main.bam", emit: bam

	script:
	def r1 = reads1
	def r2 =  file(reads2) ? reads2 : ""

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	bwa mem -a -t \$cpus ${params.reference} ${r1} ${r2} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}.main.bam -
	"""
}

process align_singles {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"

	input:
	val(sample)
	path(reads)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}.main.singles.bam", emit: bam

	script:

	"""
	cpus=\$(expr \"$task.cpus\" - 4)
	bwa mem -a -t \$cpus ${params.reference} ${reads} | samtools view -F 4 -buSh - | samtools sort -@ 2 -o ${sample}.main.singles.bam -
	"""
}

process rename_bam {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.bam", emit: bam

	script:

	"""
	mkdir -p $sample
	cp $bam ${sample}/${sample}.bam
	"""
}

process merge_and_sort {
	conda "bioconda::bwa bioconda::'samtools>=1.11'"
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(main_bam)
	path(singles_bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.bam", emit: bam

	script:
	def s_bam = file(singles_bam) ? singles_bam : ""

	if (file(singles_bam)) {
		"""
		mkdir -p $sample
		samtools merge -@ $task.cpus "${sample}/${sample}.bam" ${main_bam} ${s_bam}
		"""
	} else {
		"""
		mkdir -p $sample
		cp ${main_bam} "${sample}/${sample}.bam"
		"""
	}
}

process gffquant {
	conda params.gffquant_env
	publishDir "$output_dir", mode: params.publish_mode

	input:
	val(sample)
	path(bam)

	output:
	stdout
	val "$sample", emit: sample_id
	path "${sample}/${sample}.seqname.uniq.txt", emit: uniq_seq
	path "${sample}/${sample}.seqname.dist1.txt", emit: dist1_seq
	path "${sample}/${sample}.feature_counts.txt", emit: feat_counts
	path "${sample}/${sample}.gene_counts.txt", emit: gene_counts
	path "${sample}/${sample}.covsum.txt", emit: covsum

	script:

	"""
	mkdir -p ${sample}
	gffquant ${params.gffquant_db} ${bam} -o ${sample}/${sample} -m ${params.gffquant_mode} --ambig_mode ${params.gffquant_amode} > ${sample}/${sample}.gq.o 2> ${sample}/${sample}.gq.e
	"""
}
*/


workflow {

	/*
		Collect all fastq files from --input_dir
		Ensure your input data is gzipped and one of the following:
		- 1 single end fastq (files need to be named <sample>_R1 or <sample>_1, file endings can be .fastq.gz or .fq.gz)
		- 2 paired end fastqs (<sample>_R1/R2)
		- 2 paired end fastqs + 1 single end fastq (<sample>_R1/R2, <sample>.singles <- this is essential, otherwise they will be treated as different samples)
		If you have multiple paired end/single end data sets that are supposed to be processed together (e.g. the same library sequenced on multiple lanes)
		you need to concatenate the read files beforehand.
	*/

	reads_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq.gz,fq.gz}")
		.map { file ->
			def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
			sample = sample.replaceAll(/_R?[12]$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	/*
		Normalise fastq file names by symlinking to R1/R2.
	*/

	prepare_fastqs(reads_ch)

	/*
		Preprocess (QC) reads, then split preprocessed pairs from single-end / orphaned reads.
	*/

	qc_preprocess(prepare_fastqs.out.reads)

	qc_preprocess.out.qc_reads
		.multiMap { sample, reads ->
			single: (reads[0] != null && reads[0].name.endsWith(".fastq.gz")) ? [sample, reads[0]] : 0
			paired: (reads[1] != null && reads[1].name.endsWith(".fastq.gz") && reads[2] != null && reads[2].name.endsWith(".fastq.gz")) ? [sample, [reads[1], reads[2]]] : 0
			other: 0

		}
		.set { preprocessed_reads_ch }

	/*
		Merge paired-end reads, then split merged reads from those that failed to merge.
	*/

	qc_bbmerge(preprocessed_reads_ch.paired.filter({it != 0}))

	qc_bbmerge.out.merged_reads
		.multiMap { sample, reads ->
			merged: (reads[0] != null && reads[0].name.endsWith(".fastq.gz")) ? [sample, reads[0]] : 0
			paired: (reads[1] != null && reads[1].name.endsWith(".fastq.gz") && reads[2] != null && reads[2].name.endsWith(".fastq.gz")) ? [sample, [reads[1], reads[2]]] : 0
			other: 0
		}
		.set { merged_reads_ch }

	/*
		Redirect all unpaired reads into a common channel, then concatenate them into a single unpaired fastq file.
	*/

	single_reads_ch = merged_reads_ch.merged.filter({ it != 0 }).concat(preprocessed_reads_ch.single.filter({ it != 0 }))
		.map { sample, reads ->
			return tuple(sample.replaceAll(/.singles$/, ""), reads)
		}
        .groupTuple(sort: true)

	concat_singles(single_reads_ch)

	/*
		Route all read sets (paired and unpaired) into decontamination.
	*/

	to_decontaminate_ch = concat_singles.out.concat_reads.concat(merged_reads_ch.merged.filter({ it != 0i }))
	to_decontaminate_ch.view()

}




/* workflow {
	reads_ch = Channel
	    .fromPath(params.input_dir + "/" + params.file_pattern)
    	.map { file ->
        	def sample = file.name.replaceAll(suffix_pattern, "")
        	sample = sample.replaceAll(/_[12]$/, "")
	        return tuple(sample, file)
    	}
	    .groupTuple()
	reads_ch.view()

	aux_reads_ch = Channel
		.fromPath(params.input_dir + "/" + params.single_file_pattern)
		.map { file ->
			def sample = file.name.replaceAll(suffix_pattern, "")
			sample = sample.replaceAll(/.singles/, "")
			return tuple(sample, file)
		}
		.groupTuple()
	aux_reads_ch.view()

	if (params.r1_only) {
		qc_preprocess_singles(reads_ch)
		decontaminate_singles(qc_preprocess_singles.out.sample_id, qc_preprocess_singles.out.u_fq)
		align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
		rename_bam(align_singles.out.sample_id, align_singles.out.bam)
		gffquant(rename_bam.out.sample_id, rename_bam.out.bam)
	} else {
		qc_preprocess(reads_ch)
		decontaminate(qc_preprocess.out.sample_id, qc_preprocess.out.r1_fq, qc_preprocess.out.r2_fq)
		align(decontaminate.out.sample_id, decontaminate.out.r1_fq, decontaminate.out.r2_fq)
		if (run_singles) {
			qc_preprocess_singles(aux_reads_ch)
			merge_singles(qc_preprocess.out.sample_id, qc_preprocess.out.u_fq, qc_preprocess_singles.out.u_fq)
			decontaminate_singles(merge_singles.out.sample_id, merge_singles.out.u_fq)
			align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
			merge_and_sort(align.out.sample_id, align.out.bam, align_singles.out.bam)
		} else {
			decontaminate_singles(qc_preprocess.out.sample_id, qc_preprocess.out.u_fq)
			align_singles(decontaminate_singles.out.sample_id, decontaminate_singles.out.u_fq)
			merge_and_sort(align.out.sample_id, align.out.bam, align_singles.out.bam)
		}
		gffquant(merge_and_sort.out.sample_id, merge_and_sort.out.bam)
	}

} */
