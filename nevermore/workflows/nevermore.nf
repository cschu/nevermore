#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./prep"
include { remove_host_kraken2_individual; remove_host_kraken2 } from "../modules/decon/kraken2"
include { prepare_fastqs } from "../modules/converters/prepare_fastqs"
include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { collate_stats } from "../modules/collate"
include { nevermore_align; nevermore_prep_align } from "./align"

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


workflow nevermore_main {

	take:
		fastq_ch

	main:
		if (do_preprocessing) {
	
			//prepare_fastqs(fastq_ch)
	
			//raw_fastq_ch = prepare_fastqs.out.reads
	
			nevermore_simple_preprocessing(fastq_ch)
	
			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			if (!params.drop_orphans) {
				preprocessed_ch = preprocessed_ch.concat(nevermore_simple_preprocessing.out.orphan_reads_out)
			}
	
			if (params.remove_host) {
	
				remove_host_kraken2_individual(preprocessed_ch, params.remove_host_kraken2_db)
	
				preprocessed_ch = remove_host_kraken2_individual.out.reads
				if (!params.drop_chimeras) {
					chimera_ch = remove_host_kraken2_individual.out.chimera_orphans
						.map { sample, file ->
							def meta = [:]
							meta.id = sample.id + ".chimeras"
							meta.is_paired = false
							meta.library = sample.library
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

		if (do_preprocessing) {
			collate_ch = nevermore_simple_preprocessing.out.raw_counts
			.map { sample, file -> return file }
			.collect()
			.concat(
				nevermore_prep_align.out.read_counts
					.map { sample, file -> return file }
					.collect()
			)
		}			

		if (!params.skip_alignment) {
			nevermore_align(nevermore_prep_align.out.fastqs)
			align_ch = nevermore_align.out.alignments
	
			if (do_preprocessing) {		
				collate_ch = collate_ch	
					.concat(
						nevermore_align.out.aln_counts
							.map { sample, file -> return file }
							.collect()
					)		
			}
		}

		if (do_preprocessing) {
			collate_stats(collate_ch.collect())
		}

	emit:
		alignments = align_ch
		fastqs = nevermore_prep_align.out.fastqs
		read_counts = nevermore_prep_align.out.read_counts
}
