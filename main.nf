#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"
include { collate_stats } from "./nevermore/modules/collate"


if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir
def do_alignment = params.run_gffquant || !params.skip_alignment
def do_stream = params.gq_stream
def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


params.ignore_dirs = ""


workflow {

	fastq_input(
		Channel.fromPath(input_dir + "/*", type: "dir")
			.filter { !params.ignore_dirs.split(",").contains(it.name) },
		Channel.of(null)
	)

	fastq_ch = fastq_input.out.fastqs
	
	nevermore_main(fastq_ch)

	align_ch = Channel.empty()
	counts_ch = nevermore_main.out.readcounts

	if (!do_stream && do_alignment) {
		nevermore_align(nevermore_main.out.fastqs)
		align_ch = nevermore_align.out.alignments
		counts_ch = counts_ch.concat(
			nevermore_align.out.aln_counts
				.map { sample, file -> return file }
				.collect()
		)
	}

	if (do_preprocessing && params.run_qa) {
		collate_stats(counts_ch.collect())		
	}

	if (params.run_gffquant) {

		if (params.gq_stream) {
			gq_input_ch = nevermore_main.out.fastqs
				.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, [fastqs].flatten())
			}
			.groupTuple()
			.map { sample_id, fastqs -> return tuple(sample_id, [fastqs].flatten()) }
			gq_input_ch.view()

		} else {

			gq_input_ch = align_ch

		}

		gffquant_flow(gq_input_ch)		

	}

	if (params.run_motus) {

		motus(nevermore_main.out.fastqs, params.motus_db)

	}

}
