#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include { fastqc } from "../modules/qc/fastqc"
// include { multiqc } from "../modules/qc/multiqc"
include { bwa_mem_align } from "../modules/align/bwa"
include { merge_and_sort } from "../modules/align/helpers"
// include { merge_single_fastqs } from "../modules/converters/merge_fastqs"

def asset_dir = "${projectDir}/nevermore/assets"
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

// workflow nevermore_prep_align {

// 	take:
// 		fastq_ch
	
// 	main:
// 		/*	route all single-read files into a common channel */

// 		single_ch = fastq_ch
// 			.filter { it[0].is_paired == false }
// 			.map { sample, fastq ->
// 				def meta = sample.clone()
// 				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
// 				meta.merged = false
// 				return tuple(meta, fastq)
// 			}

// 		/*	route all paired-end read files into a common channel */

// 		paired_ch = fastq_ch
// 			.filter { it[0].is_paired == true }
// 			.map { sample, fastq ->
// 				def meta = sample.clone()
// 				meta.merged = true
// 				return tuple(meta, fastq)
// 			}

// 		/*	group all single-read files by sample and route into merge-channel */

// 		single_ch
// 			.map {
// 				sample, fastq ->
// 					sample.id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles")
// 					return tuple(sample, fastq)
// 			}
// 			.branch {
// 				single_end: it[0].library == "single"
// 				paired_end: it[0].library == "paired"
// 			}
// 		.set { single_reads_ch }

// 		def se_group_size = 3 - (params.drop_chimeras ? 1 : 0) - (params.drop_orphans ? 1 : 0)

// 		single_reads_ch.paired_end
// 			.groupTuple(sort: true, size: se_group_size, remainder: true)
// 			.branch {
// 				merge: it[1].size() > 1
// 				no_merge: true
// 			}
// 			.set { pe_singles_ch }

// 		merged_single_ch = pe_singles_ch.merge
		


// 			// .map { sample, fastq  ->
// 			// 	return tuple(sample.id, fastq)
// 			// }
// 			// .map { sample_id, files ->
// 			// 	def meta = [:]
// 			// 	meta.id = sample_id
// 			// 	meta.is_paired = false
// 			// 	meta.library = "paired"
// 			// 	meta.merged = true
// 			// 	return tuple(meta, files)
// 			// }

// 		// merged_single_ch = single_reads_ch.single_end
// 		// 	.map { sample, fastq  ->
// 		// 		return tuple(sample.id, fastq)
// 		// 	}
// 		// 	.groupTuple(sort: true)
// 		// 	.map { sample_id, files ->
// 		// 		def meta = [:]
// 		// 		meta.id = sample_id
// 		// 		meta.is_paired = false
// 		// 		meta.library = "single"
// 		// 		meta.merged = true
// 		// 		return tuple(meta, files)
// 		// 	}
// 		// 	.concat(
// 		// 		single_reads_ch.paired_end
// 		// 			.map { sample, fastq  ->
// 		// 				return tuple(sample.id, fastq)
// 		// 			}
// 		// 			.groupTuple(sort: true)
// 		// 			.map { sample_id, files ->
// 		// 				def meta = [:]
// 		// 				meta.id = sample_id
// 		// 				meta.is_paired = false
// 		// 				meta.library = "paired"
// 		// 				meta.merged = true
// 		// 				return tuple(meta, files)
// 		// 			}
// 		// 	)

// 		// merged_single_ch.view()
// 		// merged_single_ch = single_ch
// 		// 	.map { sample, fastq ->
// 		// 		return tuple(
// 		// 			sample.id.replaceAll(/.(orphans|singles|chimeras)$/, ".singles"),
// 		// 			sample.library,
// 		// 			fastq
// 		// 		)
// 		// 	}
// 		// 	.groupTuple(sort: true)
// 		// 	.map { sample_id, library, files ->
// 		// 		def meta = [:]
// 		// 		meta.id = sample_id
// 		// 		meta.is_paired = false
// 		// 		meta.library = library
// 		// 		meta.merged = true
// 		// 		return tuple(meta, files)
// 		// 	}

// 		/*	then merge single-read file groups into single files */

// 		merge_single_fastqs(merged_single_ch)

// 		/* 	take all single-read files except for the qc-survivors,
// 			concat with merged single-read files (takes care of single-end qc-survivors),
// 			concat with paired-end files,
// 			and route them into a channel for post-qc fastqc analysis
// 			(THIS IS JUST FOR STATS PURPOSES, NO WORRIES, THE SE READS ARE PROCESSED PROPERLY!)
// 		*/

// 		fastqc_in_ch = single_ch
// 			.filter { ! it[0].id.endsWith(".singles") }
// 			.map { sample, fastq ->
// 				def meta = sample.clone()
// 				meta.id = fastq.name.replaceAll(/_R1.fastq.gz$/, "")
// 				meta.merged = false
// 				return tuple(meta, fastq)
// 			}
// 			.concat(pe_singles_ch.no_merge)
// 			.concat(single_reads_ch.single_end)
// 			.concat(paired_ch)
// 			.concat(merge_single_fastqs.out.fastq)

// 		/*	perform post-qc fastqc analysis and generate multiqc report on merged single-read and paired-end sets */

// 		readcounts_ch = Channel.empty()

// 		if (params.run_qa) {

// 			fastqc(fastqc_in_ch, "qc")

// 			multiqc(
// 				fastqc.out.stats
// 					.filter { it[0].merged == true || it[0].is_paired == true }
// 					.map { sample, report -> return report }
// 					.collect(),
// 				"${asset_dir}/multiqc.config",
// 				"qc"
// 			)

// 			readcounts_ch = fastqc.out.counts
// 		}

// 		fastq_prep_ch = paired_ch
// 			.concat(single_reads_ch.single_end)
// 			.concat(pe_singles_ch.no_merge)
// 			.concat(merge_single_fastqs.out.fastq)

// 	emit:
// 		fastqs = fastq_prep_ch
// 		read_counts = readcounts_ch

// }


// align_ch = Channel.empty()
// 		if (!do_stream && do_alignment) {
// 			nevermore_align(nevermore_prep_align.out.fastqs)
// 			align_ch = nevermore_align.out.alignments
	
// 			if (do_preprocessing) {		
// 				collate_ch = collate_ch	
// 					.concat(
// 						nevermore_align.out.aln_counts
// 							.map { sample, file -> return file }
// 							.collect()
// 					)		
// 			}
// 		}

// 		if (do_preprocessing && params.run_qa) {
// 			collate_stats(collate_ch.collect())
// 		}


workflow nevermore_align {

	take:
		fastq_ch
		// readcounts_ch

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

		// if (!params.skip_preprocessing || params.run_preprocessing) && params.run_qa) {		
		// 	collate_stats(
		// 		readcounts_ch.concat(
		// 			merge_and_sort.out.flagstats.map { sample, file -> return file }
		// 		)
		// 	)

		// 	collate_ch = collate_ch	
		// 			.concat(
		// 				nevermore_align.out.aln_counts
		// 					.map { sample, file -> return file }
		// 					.collect()
		// 			)		
		// 	}
		// }

		// if (do_preprocessing && params.run_qa) {
		// 	collate_stats(collate_ch.collect())
		// }


	emit:
		alignments = merge_and_sort.out.bam
		aln_counts = merge_and_sort.out.flagstats

}

