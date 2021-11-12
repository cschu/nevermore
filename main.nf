#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fq; fq2bam; prepare_fastqs } from "./modules/nevermore/convert"
include { nevermore_simple_preprocessing } from "./workflows/nevermore/nevermore"
include { classify_sample } from "./modules/nevermore/functions"
include { remove_host_kraken2 } from "./modules/nevermore/decon/kraken2"
include { flagstats; count_reads } from "./modules/nevermore/stats"

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)


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

		raw_fastq_ch = prepare_fastqs.out.reads.concat(bam2fq.out.reads)

		nevermore_simple_preprocessing(raw_fastq_ch)


		if (params.remove_host) {

			remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

			preprocessed_ch = remove_host_kraken2.out.reads

		} else {

			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out

		}

	} else {

		preprocessed_ch = fastq_ch

	}


	if (run_count_reads || run_bam_analysis) {

		fq2bam(preprocessed_ch)

		if (run_count_reads) {

	        flagstats(fq2bam.out.reads)

    	    count_reads(flagstats.out.flagstats)

		}

		if (run_bam_analysis) {

			bam_analysis(fq2bam.out.reads)

		}

    }


	if (run_fastq_analysis) {

		fastq_analysis(preprocessed_ch)

	}
}
