#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { bwa_mem_align } from "../modules/align/bwa"
include { merge_and_sort } from "../modules/align/helpers"
include { merge_single_fastqs } from "../modules/converters/merge_fastqs"

def asset_dir = "${projectDir}/nevermore/assets"


workflow nevermore_prep_align {

	take:
		fastq_ch
	
	main:
		/*	route all single-read files into a common channel */

		single_ch = fastq_ch
			.filter { it[0].is_paired == false }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}

		/*	route all paired-end read files into a common channel */

		paired_ch = fastq_ch
			.filter { it[0].is_paired == true }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = sample.id
				meta.is_paired = true
				meta.merged = true
				return tuple(meta, fastq)
			}

		/*	group all single-read files by sample and route into merge-channel */

		merged_single_ch = single_ch
			.map { sample, fastq ->
				return tuple(
					sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles"),
					fastq
				)
			}
			.groupTuple(sort: true)
			.map { sample_id, files ->
				def meta = [:]
				meta.id = sample_id
				meta.is_paired = false
				meta.merged = true
				return tuple(meta, files)
			}

		/*	then merge single-read file groups into single files */

		merge_single_fastqs(merged_single_ch)

		/* 	take all single-read files except for the qc-survivors,
			concat with merged single-read files (takes care of single-end qc-survivors),
			concat with paired-end files,
			and route them into a channel for post-qc fastqc analysis
			(THIS IS JUST FOR STATS PURPOSES, NO WORRIES, THE SE READS ARE PROCESSED PROPERLY!)
		*/

		fastqc_in_ch = single_ch
			.filter { ! it[0].id.endsWith(".singles") }
			.map { sample, fastq ->
				def meta = [:]
				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
				meta.is_paired = false
				meta.merged = false
				return tuple(meta, fastq)
			}
			.concat(merge_single_fastqs.out.fastq)
			.concat(paired_ch)

		/*	perform post-qc fastqc analysis and generate multiqc report on merged single-read and paired-end sets */

		fastqc(fastqc_in_ch, "qc")

		multiqc(
			fastqc.out.stats
				.filter { it[0].merged == true || it[0].is_paired == true }
				.map { sample, report -> return report }
				.collect(),
			"${asset_dir}/multiqc.config",
			"qc"
		)

		fastq_prep_ch = paired_ch.concat(merge_single_fastqs.out.fastq)

	emit:
		fastqs = fastq_prep_ch
		read_counts = fastqc.out.counts
		

}


workflow nevermore_align {

	take:

		fastq_ch

	main:
		
		/*	align merged single-read and paired-end sets against reference */

		bwa_mem_align(
			fastq_ch,
			params.reference,
			true
		)

		/*	merge paired-end and single-read alignments into single per-sample bamfiles */

		aligned_ch = bwa_mem_align.out.bam
			.map { sample, bam ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, bam)
			}
			.groupTuple(sort: true)

		merge_and_sort(aligned_ch, true)


	emit:

		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}

