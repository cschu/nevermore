include { stream_gffquant as gffquant; collate_feature_counts } from "../modules/profilers/gffquant"

params.gq_collate_columns = "uniq_scaled,combined_scaled"


workflow gffquant_flow {

	take:

		input_ch

	main:
		gffquant(input_ch)
		feature_count_ch = (params.gq_panda) ? gffquant.out.profiles : gffquant.out.results
		counts = gffquant.out.results		

		feature_count_ch = feature_count_ch
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith("Counter.txt.gz") }
			.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
			.filter { params.collate_gene_counts || !it.name.endsWith("gene_counts.pd.txt")}
			.filter { !it.name.endsWith("gene_ids.txt.gz") }
			.map { file -> 
				def category = file.name
					// .replaceAll(/\.txt(\.gz)?$/, "")
					.replaceAll(/\.pd\.txt$/, "")
					.replaceAll(/\.txt\.gz$/, "")
					.replaceAll(/.+\./, "")
				return [ category, file ]
			}
			.groupTuple(sort: true)
			.combine(
				Channel.from(params.gq_collate_columns.split(","))
			)

		collate_feature_counts(
			feature_count_ch,
			(params.future_features ? ((params.gq_panda) ? ".pd.txt" : ".txt.gz") : "")
		)

	emit:

		counts
		collated = collate_feature_counts.out.collated

}
