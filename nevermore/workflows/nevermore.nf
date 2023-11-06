#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./prep"
include { remove_host_kraken2_individual; remove_host_kraken2 } from "../modules/decon/kraken2"
include { sortmerna } from "../modules/decon/sortmerna"
include { prepare_fastqs } from "../modules/converters/prepare_fastqs"
include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { collate_stats } from "../modules/collate"
include { nevermore_align; nevermore_prep_align } from "./align"


workflow nevermore_main {

	take:
		fastq_ch
		

	main:
		if (params.qc.run) {
	
			nevermore_simple_preprocessing(fastq_ch)
	
			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			if (params.qc.keep_orphans) {
				preprocessed_ch = preprocessed_ch.concat(nevermore_simple_preprocessing.out.orphan_reads_out)
			}

			if (params.run_sortmerna) {

				preprocessed_ch
					.branch {
						metaT: it[0].library_type == "metaT"
						metaG: true
					}
					.set { for_sortmerna_ch }

				sortmerna(for_sortmerna_ch.metaT, params.sortmerna_db)
				preprocessed_ch = for_sortmerna_ch.metaG
					.concat(sortmerna.out.fastqs)

			}
	
			if (params.decon.run) {
	
				remove_host_kraken2_individual(preprocessed_ch, params.decon.kraken2.db)
	
				preprocessed_ch = remove_host_kraken2_individual.out.reads
				if (params.decon.keep_chimeras) {
					chimera_ch = remove_host_kraken2_individual.out.chimera_orphans
						.map { sample, file ->
							def meta = sample.clone()
							meta.id = sample.id + ".chimeras"
							meta.is_paired = false
							return tuple(meta, file)
						}
					preprocessed_ch = preprocessed_ch.concat(chimera_ch)
				}
	
			}
	
		} else {
	
			preprocessed_ch = fastq_ch
	
		}
	
		nevermore_prep_align(preprocessed_ch)
		align_ch = Channel.empty()
		collate_ch = Channel.empty()

		if (params.qc.run) {
			collate_ch = nevermore_simple_preprocessing.out.raw_counts
			.map { sample, file -> return file }
			.collect()
			.concat(
				nevermore_prep_align.out.read_counts
					.map { sample, file -> return file }
					.collect()
			)
		}

		if (params.align.run || (params.profilers.gffquant.run && !params.profilers.gffquant.stream)) {
			nevermore_align(nevermore_prep_align.out.fastqs)
			align_ch = nevermore_align.out.alignments
	
			if (params.qc.run) {		
				collate_ch = collate_ch	
					.concat(
						nevermore_align.out.aln_counts
							.map { sample, file -> return file }
							.collect()
					)		
			}
		}

		if (params.qc.run && params.qc.generate_reports) {
			collate_stats(collate_ch.collect())
		}

	emit:
		alignments = align_ch
		fastqs = nevermore_prep_align.out.fastqs
		read_counts = nevermore_prep_align.out.read_counts
}
