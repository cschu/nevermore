process run_gffquant {
	label "gffquant"

	input:
	tuple val(sample), path(alignments)
	path(gq_db)

	output:
	tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results
	tuple val(sample), path("logs/${sample}.log")

	script:
	def gq_output = "-o profiles/${sample}/${sample}"

	def gq_params = "-m ${params.gq_mode} --ambig_mode ${params.gq_ambig_mode}"
	gq_params += (params.gq_strand_specific) ? " --strand_specific" : ""
	gq_params += (params.gq_unmarked_orphans) ? " --unmarked_orphans" : ""
	gq_params += (params.gq_calc_coverage) ? " --calc_coverage" : ""
	gq_params += (params.gq_min_seqlen) ? (" --min_seqlen " + params.gq_min_seqlen) : ""
	gq_params += (params.gq_min_identity) ? (" --min_identity " + params.gq_min_identity) : ""
	gq_params += (params.bam_input_pattern || !params.large_reference) ? (" --format bam") : " --format sam"

	def gq_cmd = "gffquant ${gq_output} ${gq_params} gq_db.sqlite3"

	def mk_aln_sam = ""
	if (params.bam_input_pattern) {

		if (params.do_name_sort) {
			gq_cmd = "samtools collate -@ ${task.cpus} -O ${alignments} tmp/collated_bam | ${gq_cmd} -"
		} else {
			gq_cmd = "${gq_cmd} ${alignments}"
		}

	} else if (params.large_reference) {

		mk_aln_sam += "echo 'Making alignment stream...'\n"
		if (alignments instanceof Collection && alignments.size() >= 2) {
			mk_aln_sam += "cat ${sample}.sam > tmp/alignments.sam \n"
			mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> tmp/alignments.sam"
		} else {
			mk_aln_sam += "ln -s ${alignments[0]} tmp/alignments.sam"
		}
		gq_cmd = "cat tmp/alignments.sam | ${gq_cmd} -"

	} else {

		gq_cmd = "${gq_cmd} ${alignments}"

	}

	"""
	set -e -o pipefail
	mkdir -p logs/ tmp/ profiles/
	echo 'Copying database...'
	cp -v ${gq_db} gq_db.sqlite3
	${mk_aln_sam}
	${gq_cmd} &> logs/${sample}.log
	rm -rfv gq_db.sqlite3* tmp/
	"""
}

// process collate_feature_counts {

// 	input:
// 	tuple val(sample), path(count_tables)

// 	output:
// 	path("profiles/collated/*.txt.gz"), emit: collated, optional: true

// 	script:
// 	"""
// 	mkdir -p profiles/collated/
// 	collate_counts . -o profiles/collated/collated -c uniq_scaled
// 	collate_counts . -o profiles/collated/collated -c combined_scaled
// 	"""
// }

params.gq_collate_columns = "uniq_scaled,combined_scaled"

process collate_feature_counts {
	// publishDir "${params.output_dir}", mode: params.publish_mode

	input:
	tuple val(sample), path(count_tables), val(column)

	output:
	path("collated/*.txt.gz"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/

	collate_counts . -o collated/collated -c ${column}
	"""
}

	// for column in \$(echo ${params.gq_collate_columns} | sed 's/,/ /g'); do
	// 	collate_counts . -o collated/collated -c \$column
	// done