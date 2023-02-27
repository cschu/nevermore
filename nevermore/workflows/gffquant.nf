include { run_gffquant; collate_feature_counts; } from "../modules/profilers/gffquant"


workflow gffquant_flow {

	take:

		bam_ch

	main:

		run_gffquant(bam_ch, params.gffquant_db)

		feature_count_ch = run_gffquant.out.results
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith("Counter.txt.gz") }
			.map { file -> 
				def category = file.name
					.replaceAll(/\.txt\.gz$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort: true)

		collate_feature_counts(feature_count_ch)

	emit:

		counts = run_gffquant.out.results
		collated = collate_feature_counts.out.collated

}
