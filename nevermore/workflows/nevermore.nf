#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_simple_preprocessing } from "./prep"
include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"
include { collate_stats } from "../modules/stats"
include { nevermore_align } from "./align"
include { nevermore_pack_reads } from "./pack"
include { nevermore_qa } from "./qa"
include { nevermore_decon } from "./decon"


params.run_preprocessing = params.run_qc
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream


workflow nevermore_main {

	take:
		fastq_ch

	main:
		if (do_preprocessing) {
	
			nevermore_simple_preprocessing(fastq_ch)
	
			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			if (!params.drop_orphans) {
				preprocessed_ch = preprocessed_ch.mix(nevermore_simple_preprocessing.out.orphan_reads_out)
			}

			nevermore_decon(preprocessed_ch)
			preprocessed_ch = nevermore_decon.out.reads

		} else {
	
			preprocessed_ch = fastq_ch
	
		}
	
		nevermore_pack_reads(preprocessed_ch)

		collate_ch = Channel.empty()
		if (params.run_qa) {

			raw_counts_ch = (do_preprocessing) ? nevermore_simple_preprocessing.out.raw_counts : Channel.empty()

			nevermore_qa(
				nevermore_pack_reads.out.qa_fastqs,
				raw_counts_ch
			)

			collate_ch = nevermore_qa.out.readcounts_ch

		}

		collate_stats(collate_ch.collect())


		preprocessed_fastq_ch = nevermore_pack_reads.out.fastqs
			.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/\.singles$/, "")
				return [ sample_id, sample.is_paired, fastqs ]  //tuple(sample_id, fastqs)
			}
			.groupTuple(size: ((params.single_end_libraries) ? 1 : 2), remainder: true)
			.map { sample_id, pair_info, fastqs ->
				def meta = [:]
				meta.id = sample_id
				// meta.is_paired = pair_info.contains(true)
				return tuple(meta, [fastqs].flatten())
			}

	emit:
		// fastqs = nevermore_pack_reads.out.fastqs
		fastqs = preprocessed_fastq_ch
		readcounts = collate_ch

}
