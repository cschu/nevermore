params.profilers.gffquant.aligner = "bwa"
params.profilers.gffquant.collate_columns = "uniq_scaled,combined_scaled"


def compile_param_string(sample_id, cpus, bam_input) {
	def param_str = "-m ${params.profilers.gffquant.mode} --ambig_mode ${params.profilers.gffquant.ambig_mode}"
	param_str += (params.profilers.gffquant.strand_specific) ? " --strand_specific" : ""
	param_str += (params.profilers.gffquant.min_seqlen) ? (" --min_seqlen " + params.profilers.gffquant.min_seqlen) : ""
	param_str += (params.profilers.gffquant.min_identity) ? (" --min_identity " + params.profilers.gffquant.min_identity) : ""
	param_str += (params.profilers.gffquant.restrict_metrics) ? " --restrict_metrics ${params.profilers.gffquant.restrict_metrics}" : ""
	param_str += " -t ${cpus}"

	if (params.profilers.gffquant.gq_mode == "domain") {
		param_str += " --db_separator , --db_coordinates hmmer"
	}

	if (bam_input) {
		param_str += (params.profilers.gffquant.unmarked_orphans) ? " --unmarked_orphans" : ""
		param_str += (params.input.bam_input_pattern || !params.input.large_reference) ? (" --format bam") : " --format sam" // not sure if that still works with recent gffquant versions?
	} else {
		param_str += (params.profilers.gffquant.keep_alignments) ? " --keep_alignment_file ${sample_id}.sam" : ""
	}

	return param_str
}


process stream_gffquant {
	label "gffquant"
	tag "gffquant.${sample}"

	input:
		tuple val(sample), path(fastqs)
		path(gq_db)
		path(reference)
	output:
		tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results
		tuple val(sample), path("logs/${sample}.log")

	script:
			def gq_output = "-o profiles/${sample}/${sample}"

			def gq_params = compile_param_string(sample, task.cpus, false)

			def input_files = ""
			// we cannot auto-detect SE vs. PE-orphan! --> i think this can be read from the sample object TODO!
			if (params.profilers.gffquant.single_end_library) {
				input_files += "--fastq-singles ${fastqs}"
			} else {
				r1_files = fastqs.findAll( { it.name.endsWith("_R1.fastq.gz") && !it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )
				r2_files = fastqs.findAll( { it.name.endsWith("_R2.fastq.gz") } )
				orphans = fastqs.findAll( { it.name.matches("(.*)(singles|orphans|chimeras)(.*)") } )

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
	
			def gq_cmd = "gffquant ${gq_output} ${gq_params} --db GQ_DATABASE --reference \$(readlink ${reference}) --aligner ${params.profilers.gffquant.aligner} ${input_files}"

			"""
			set -e -o pipefail
			mkdir -p logs/ tmp/ profiles/
			echo 'Copying database...'
			cp -v ${gq_db} GQ_DATABASE
			${gq_cmd} &> logs/${sample}.log
			rm -rfv GQ_DATABASE* tmp/
			"""

}

process run_gffquant {
	label "gffquant"
	tag "gffquant.${sample}"

	input:
	tuple val(sample), path(alignments)
	path(gq_db)

	output:
	tuple val(sample), path("profiles/${sample}/*.txt.gz"), emit: results
	tuple val(sample), path("logs/${sample}.log")

	script:
	def gq_output = "-o profiles/${sample}/${sample}"

	def gq_params = compile_param_string(sample, task.cpus, true)

	def gq_cmd = "gffquant ${gq_output} ${gq_params} GQ_DATABASE"

	def mk_aln_sam = ""
	if (params.input.bam_input_pattern && params.profilers.gffquant.do_name_sort) {

		gq_cmd = "samtools collate -@ ${task.cpus} -O ${alignments} tmp/collated_bam | ${gq_cmd} --bam -"

	} else if (params.input.large_reference) {

		mk_aln_sam += "echo 'Making alignment stream...'\n"
		if (alignments instanceof Collection && alignments.size() >= 2) {
			mk_aln_sam += "cat ${sample}.sam > tmp/alignments.sam \n"
			mk_aln_sam += "grep -v '^@' ${sample}.singles.sam >> tmp/alignments.sam"
		} else {
			mk_aln_sam += "ln -s ${alignments[0]} tmp/alignments.sam"
		}
		gq_cmd = "cat tmp/alignments.sam | ${gq_cmd} --sam -"

	} else {

		gq_cmd = "${gq_cmd} --bam ${alignments}"

	}

	"""
	set -e -o pipefail
	mkdir -p logs/ tmp/ profiles/
	echo 'Copying database...'
	cp -v ${gq_db} GQ_DATABASE
	${mk_aln_sam}
	${gq_cmd} &> logs/${sample}.log
	rm -rfv GQ_DATABASE* tmp/
	"""
}



process collate_feature_counts {

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
