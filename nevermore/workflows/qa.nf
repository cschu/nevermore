include { fastqc } from "../modules/qc/fastqc"
include { multiqc } from "../modules/qc/multiqc"


workflow nevermore_qa {

	take:
		fastq_ch
		counts_ch

	main:
		fastqc(fastqc, "qc")

		multiqc(
			fastqc.out.stats
				.filter { it[0].merged == true || it[0].is_paired == true }
				.map { sample, report -> return report }
				.collect(),
			"${asset_dir}/multiqc.config",
			"qc"
		)

		readcounts_ch = counts_ch
			.concat(fastqc.out.counts)
			.map { sample, file -> return file }	

	emit:
		readcounts_ch

}
