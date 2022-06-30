#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./nevermore/workflows/nevermore"
include { classify_sample } from "./nevermore/modules/functions"
include { remove_host_kraken2_individual; remove_host_kraken2 } from "./nevermore/modules/decon/kraken2"
include { prepare_fastqs } from "./nevermore/modules/converters/prepare_fastqs"

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


process bwa_mem_align {
	label 'align'

	input:
	tuple val(sample), path(reads)
	path(reference)

	output:
	tuple val(sample), path("${sample.id}.bam"), emit: bam

	script:
	def align_cpus = 4 // figure out the groovy division garbage later (task.cpus >= 8) ?
	def sort_cpus = 4
	def reads2 = (sample.is_paired) ? "${sample.id}_R2.fastq.gz" : ""

	"""
	bwa mem -a -t ${align_cpus} \$(readlink ${reference}) ${sample.id}_R1.fastq.gz ${reads2} | samtools view -F 4 -buSh - | samtools sort -@ ${sort_cpus} -o ${sample.id}.bam
	"""

}

process merge_and_sort {
	label 'samtools'
	publishDir params.output_dir, mode: params.publish_mode

	input:
	tuple val(sample), path(bamfiles)

	output:
	tuple val(sample), path("bam/${sample}.bam"), emit: bam

	script:
	// need a better detection for this
	if (bamfiles instanceof Collection && bamfiles.size() >= 2) {
		"""
		mkdir -p bam/
		samtools merge -@ $task.cpus bam/${sample}.bam ${bamfiles}
		"""
	} else {
		// i don't like this solution
		"""
		mkdir -p bam
		ln -s ${bamfiles[0]} bam/${sample}.bam
		"""
	}
}

process gffquant {
	publishDir params.output_dir, mode: params.publish_mode

	input:
	tuple val(sample), path(bam)
	path(annotation_db)

	output:
	path "${sample}/${sample}.seqname.uniq.txt", emit: uniq_seq
	path "${sample}/${sample}.seqname.dist1.txt", emit: dist1_seq
	path "${sample}/${sample}.feature_counts.txt", emit: feat_counts
	path "${sample}/${sample}.gene_counts.txt", emit: gene_counts
	path "${sample}/${sample}.covsum.txt", emit: covsum

	script:

	"""
	mkdir -p ${sample}
	gffquant ${annotation_db} ${bam} -o ${sample}/${sample} -m ${params.gffquant_mode} --ambig_mode ${params.gffquant_amode}
	"""
}



/*
workflow {
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

}
*/


workflow {

	fastq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
		.map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)
        .map { classify_sample(it[0], it[1]) }


	if (do_preprocessing) {

		prepare_fastqs(fastq_ch)

		raw_fastq_ch = prepare_fastqs.out.reads

		nevermore_simple_preprocessing(raw_fastq_ch)

		preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
		if (!params.drop_orphans) {
			preprocessed_ch = preprocessed_ch.concat(nevermore_simple_preprocessing.out.orphan_reads_out)
		}

		if (params.remove_host) {

			// remove_host_kraken2(preprocessed_ch, params.remove_host_kraken2_db)
			remove_host_kraken2_individual(preprocessed_ch, params.remove_host_kraken2_db)

			// preprocessed_ch = remove_host_kraken2.out.reads
			preprocessed_ch = remove_host_kraken2_individual.out.reads
			if (!params.drop_chimeras) {
				chimera_ch = remove_host_kraken2_individual.out.chimera_orphans
					.map { sample, file ->
						def meta = [:]
						meta.is_paired = false
						meta.id = sample.id + ".chimeras"
						return tuple(meta, file)
					}
				preprocessed_ch = preprocessed_ch.concat(chimera_ch)
			}

		}

	} else {

		preprocessed_ch = fastq_ch

	}

	bwa_mem_align(preprocessed_ch, params.reference)

	aligned_ch = bwa_mem_align.out.bam
		.map { sample, bam ->
			sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
			return tuple(sample_id, bam)
		}
		.groupTuple(sort: true)

	aligned_ch.view()

	merge_and_sort(aligned_ch)
	//gffquant(merge_and_sort.out.bam, params.gffquant_db)
}
