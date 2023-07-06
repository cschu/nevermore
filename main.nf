#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

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

params.ignore_dirs = ""



workflow {

	fastq_input(
		Channel.fromPath(input_dir + "/*", type: "dir")
			.filter { !params.ignore_dirs.split(",").contains(it.name) }
	)

	fastq_ch = fastq_input.out.fastqs
	
	nevermore_main(fastq_ch)

	if (params.run_gffquant) {

		if (params.gq_stream) {
			gq_input_ch = nevermore_main.out.fastqs
				.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, [fastqs].flatten())
			}
			.groupTuple()
			.map { sample_id, fastqs -> return tuple(sample_id, fastqs.flatten()) }
			gq_input_ch.view()
			//.groupTuple(sort: true)

		} else {

			gq_input_ch = nevermore_main.out.alignments

		}

		// gq_input_ch = ((params.gq_stream) ? nevermore_main.out.fastqs : never_main.out.alignments)
		// 	.map { sample, files -> return tuple(sample.id, files) }
		gffquant_flow(gq_input_ch)		

	}

	if (params.run_motus) {

		motus(nevermore_main.out.fastqs, params.motus_db)

	}

}
