params.gq_aligner = "bwa"
params.gq_min_seqlen = params.min_alignment_length
params.gq_single_end_library = params.single_end_libraries
params.gq_min_identity = params.min_identity
params.gq_mode = "genes"
params.gq_ambig_mode = "1overN"
params.gq_restrict_metrics = "raw,lnorm,scaled,rpkm"


process stream_gffquant {
	container "ghcr.io/cschu/gff_quantifier:v0.20.0"
	tag "gffquant.${sample}"
	label "gffquant"
	label "large"

	input:
		tuple val(sample), path(fastqs), path(gq_db)

	output:
		tuple val(sample), path("${sample}/*.{txt.gz,pd.txt}"), emit: results
		tuple val(sample), path("${sample}/*.{txt.gz,pd.txt}"), emit: profiles
		tuple val(sample), path("logs/${sample}.log")
		tuple val(sample), path("alignments/${sample}/${sample}*.sam"), emit: alignments, optional: true
		tuple val(sample), path("${sample}/${sample}.gene_ids.txt.gz"), emit: gene_ids

	script:
			def gq_output = "-o ${sample}/${sample}"

			def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
			gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
			gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
			gq_params += (params.gq_gene_group_db) ? " --gene_group_db" : ""
			// LEGACY PARAMETERS, partially not implemented in newer gffquant
			gq_params += (params.gq_restrict_metrics) ? " --restrict_metrics ${params.gq_restrict_metrics}" : ""
			def mkdir_alignments = (params.keep_alignment_file != null && params.keep_alignment_file != false) ? "mkdir -p alignments/${sample}/" : ""

			gq_params += " -t ${task.cpus}"

			if (params.gq_mode == "domain") {
				gq_params += " --db_separator , --db_coordinates hmmer"
			}

			def input_files = ""
			r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
			r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
			orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

			if (params.gq_single_end_library || (r1_files.size() + r2_files.size() + orphans.size() == 1)) {

				input_files += "--fastq-singles ${fastqs}"

			} else {

				if (r1_files.size() != 0) {
					input_files += "--fastq-r1 ${r1_files.join(' ')}"
				}
				if (r2_files.size() != 0) {
					input_files += " --fastq-r2 ${r2_files.join(' ')}"
				}
				if (orphans.size() != 0) {
					input_files += " --fastq-orphans ${orphans.join(' ')}"
				}

			}
	

			def db_command = (params.copy_database) 
				? "cp -v \$(dirname \$(readlink ${gq_db}))/*.sqlite3 GQ_DATABASE"
				: "GQ_DATABASE=\$(dirname \$(readlink ${gq_db}))/*.sqlite3"
			
			def db_call = (params.copy_database) ? "--db GQ_DATABASE" : "--db \$GQ_DATABASE"
			def db_clean = (params.copy_database) ? "rm -fv GQ_DATABASE" : ""
			
			// def gq_cmd = "gffquant ${gq_output} ${gq_params} --db \$GQ_DATABASE --aligner ${params.gq_aligner} ${input_files}"
			def gq_cmd = "gffquant ${gq_output} ${gq_params} ${db_call} --aligner ${params.gq_aligner} ${input_files}"
			// GQ_DATABASE=\$(dirname \$(readlink ${gq_db}))/*sqlite3
			"""
			set -e -o pipefail
			mkdir -p logs/ tmp/
			${mkdir_alignments}
			${db_command}

			${gq_cmd} --reference \$(readlink ${gq_db}) | tee logs/${sample}.log
			gzip -dc ${sample}/${sample}.gene_counts.txt.gz | cut -f 1 | gzip -c - > ${sample}/${sample}.gene_ids.txt.gz
			rm -rfv tmp/
			${db_clean}
			"""

}

params.gq_collate_columns = "uniq_scaled,combined_scaled"

process collate_feature_counts {
	publishDir params.output_dir, mode: "copy"
	label "collate_profiles"
	label "gffquant"

	input:
	tuple val(sample), path(count_tables), val(column)
	val(suffix)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	def suffix_param = (suffix != "") ? "--suffix ${suffix}" : ""
	"""
	mkdir -p collated/

	collate_counts . -o collated/collated -c ${column} ${suffix_param}
	"""
}
