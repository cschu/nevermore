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

	if (!params.skip_profiling) {

		gffquant_flow(nevermore_main.out.alignments)

	}

}
