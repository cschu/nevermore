include { run_samestr_convert; run_samestr_merge; run_samestr_filter; run_samestr_stats; run_samestr_compare; run_samestr_summarize } from "../modules/profilers/samestr"

params.samestr_marker_db = "/scratch/schudoma/databases/samestr/mpa_vOct22_CHOCOPhlAnSGB_202212/marker_db/"

workflow samestr {
    take:
        samestr_convert_ch

    main:
        run_samestr_convert(
			samestr_convert_ch,
			params.samestr_marker_db
		)

        Channel
            run_samestr_convert.out.sstr_npy
            .flatten()
            .map { file ->
                    def species = file.name.replaceAll(/[.].*/, "")
                    return tuple(species, file)
            }
            .groupTuple(sort: true)
            .set { grouped_npy_ch }

		run_samestr_merge(grouped_npy_ch)
		run_samestr_filter(
			run_samestr_merge.out.sstr_npy
			params.samestr_marker_db
		)
		run_samestr_stats(run_samestr_filter.out.sstr_npy)
		run_samestr_compare(run_samestr_filter.out.sstr_npy)

		// symlink all sstr_compare/mp_profiles
        run_samestr_summarize(
			run_samestr_compare.out.sstr_compare.collect(),
			samestr_convert_ch
				.map { sam, profile -> return profile }
				.collect()
		)
 
    emit:
        run_samestr_summarize.out.sstr_summarize
        run_samestr_stats.out.sstr_stats
        run_samestr_filter.out.sstr_npy
}