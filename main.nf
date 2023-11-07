#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"
include { gffquant_flow } from "./nevermore/workflows/gffquant"
include { fastq_input } from "./nevermore/workflows/input"

include { get_single_input_dir } from "./nevermore/modules/functions/validation"

def input_dir = get_single_input_dir()

params.input.ignore_dirs = null


workflow {

	fastq_input(
		Channel.fromPath(input_dir + "/*", type: "dir")
			.filter { params.input.ignore_dirs and !params.input.ignore_dirs.split(",").contains(it.name) }
	)

	fastq_ch = fastq_input.out.fastqs
	
	nevermore_main(fastq_ch)

	if (params.profilers.gffquant.run) {

		gffquant_flow((params.profilers.gffquant.stream) ? nevermore_main.out.fastqs : nevermore_main.out.alignments)

	}

	if (params.profilers.motus.run) {

		motus(nevermore_main.out.fastqs, params.profilers.motus.db)

	}

}
