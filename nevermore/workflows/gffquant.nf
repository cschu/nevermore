include { stream_gffquant; run_gffquant; collate_feature_counts } from "../modules/profilers/gffquant"

params.profilers.gffquant.collate_columns = "uniq_scaled,combined_scaled"
params.profilers.gffquant.collate_gene_counts = true

workflow gffquant_flow {

	take:

		input_ch

	main:

		if (params.profilers.gffquant.stream) {
			gq_input_ch = input_ch
				.map { sample, fastqs ->
				sample_id = sample.id.replaceAll(/.(orphans|singles|chimeras)$/, "")
				return tuple(sample_id, [fastqs].flatten())
			}
			.groupTuple()
			.map { sample_id, fastqs -> return tuple(sample_id, fastqs.flatten()) }
			
			stream_gffquant(gq_input_ch, params.profilers.gffquant.db, params.profilers.gffquant.reference_index)
			feature_count_ch = stream_gffquant.out.results
			counts = stream_gffquant.out.results

		} else {

			run_gffquant(input_ch, params.profilers.gffquant.db)
			feature_count_ch = run_gffquant.out.results
			counts = run_gffquant.out.results

		}

		feature_count_ch = feature_count_ch
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith("Counter.txt.gz") }
			.filter { params.profilers.gffquant.collate_gene_counts || !it.name.endsWith("gene_counts.txt.gz") }
			.map { file -> 
				def category = file.name
					.replaceAll(/\.txt\.gz$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort: true)
			.combine(
				Channel.from(params.profilers.gffquant.collate_columns.split(","))
			)

		collate_feature_counts(feature_count_ch)

	emit:

		counts
		collated = collate_feature_counts.out.collated

}
