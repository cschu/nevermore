include { run_gffquant; collate_feature_counts; } from "../modules/profilers/gffquant"


workflow gffquant_flow {

	take:

		bam_ch

	main:

		run_gffquant(bam_ch, params.gffquant_db)

		feature_count_ch = run_gffquant.out.results
			.map { sample, files -> return files }
			.flatten()
			.filter { !it.name.endsWith(".aln_stats.txt.gz") }
			.filter { !it.name.endsWith("Counter.txt.gz") }
			.map { file -> 
				def category = file.name
					.replaceAll(/\.txt\.gz$/, "")
					.replaceAll(/.+\./, "")
				return tuple(category, file)
			}
			.groupTuple(sort: true)

		// M0x20MCx1884.aln_stats.txt.gz            M0x20MCx1884.BRITE.txt.gz                    M0x20MCx1884.EC_number.txt.gz    M0x20MCx1884.Gene_Ontology_terms.txt.gz  M0x20MCx1884.KEGG_Pathway.txt.gz   M0x20MCx1884.KEGG_TC.txt.gz
		// M0x20MCx1884.AmbiguousSeqCounter.txt.gz  M0x20MCx1884.CAZy.txt.gz                     M0x20MCx1884.eggNOG_OGs.txt.gz   M0x20MCx1884.KEGG_ko.txt.gz              M0x20MCx1884.KEGG_rclass.txt.gz    M0x20MCx1884.UniqueSeqCounter.txt.gz
		// M0x20MCx1884.BiGG_Reaction.txt.gz        M0x20MCx1884.COG_Functional_Category.txt.gz  M0x20MCx1884.gene_counts.txt.gz  M0x20MCx1884.KEGG_Module.txt.gz          M0x20MCx1884.KEGG_Reaction.txt.gz

		// feature_count_ch = run_gffquant.out.results //.collect()
		// 	.map { sample, files -> return files }
		// 	.flatten()
		// 	.filter { !it.name.endsWith("gene_counts.txt") }
		// 	.filter { !it.name.endsWith("seqname.uniq.txt") }
		// 	.filter { !it.name.endsWith("seqname.dist1.txt") }
		// 	.map { file -> 
		// 		def category = file.name.replaceAll(/\.txt$/, "")
		// 			.replaceAll(/.+\./, "")
		// 		return tuple(category, file)
		// 	}
		// 	.groupTuple(sort:true)

		collate_feature_counts(feature_count_ch)

	emit:

		counts = run_gffquant.out.results
		collated = collate_feature_counts.out.collated

}
