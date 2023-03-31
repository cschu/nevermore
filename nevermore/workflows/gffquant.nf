include { run_gffquant; collate_feature_counts; } from "../modules/profilers/gffquant"

params.gq_collate_columns = "uniq_scaled,combined_scaled"


workflow gffquant_flow {

	take:

		bam_ch

	main:

		run_gffquant(bam_ch, params.gffquant_db)

		feature_count_ch = run_gffquant.out.results
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith("Counter.txt.gz") }
			.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
			.map { file -> 
				def category = file.name
					.replaceAll(/\.txt\.gz$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort: true)
			.combine(
				Channel.from(params.gq_collate_columns.split(","))
			)

		collate_feature_counts(feature_count_ch)

	emit:

		counts = run_gffquant.out.results
		collated = collate_feature_counts.out.collated

}
