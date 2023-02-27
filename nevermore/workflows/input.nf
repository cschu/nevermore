nextflow.enable.dsl=2

include { classify_sample } from "../modules/functions"


if (!params.bam_input_pattern) {
	params.bam_input_pattern = "**.bam"
}

def bam_suffix_pattern = params.bam_input_pattern.replaceAll(/\*/, "")

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir


process transfer_fastqs {
	input:
		path(fastqs)
	output:
		path("fastq/*/**"), emit: fastqs
	script:
		"""
		transfer.py -i . -o fastq/
		"""
}


process transfer_bams {
	input:
		path(bamfiles)
	output:
		path("bam/*.bam"), emit: bamfiles
	script:
		"""
		mkdir -p bam/
		find . -maxdepth 1 -type l -name '*.bam' | xargs -I {} readlink {} | xargs -I {} rsync -vP {} bam/
		"""
}


process prepare_fastqs {

	input:
		path(files)
		val(remote_input)
	output:
		path("fastq/*/*.fastq.{gz,bz2}"), emit: fastqs

  script:
		def remote_option = (remote_input) ? "--remote-input" : ""
		def remove_suffix = (params.suffix_pattern) ? "--remove-suffix ${params.suffix_pattern}" : ""
		def input_dir_prefix = (params.input_dir) ? params.input_dir : params.remote_input_dir

		def custom_suffixes = (params.custom_fastq_file_suffixes) ? "--valid-fastq-suffixes ${params.custom_fastq_file_suffixes}" : ""
		
		"""
		prepare_fastqs.py -i . -o fastq/ -p ${input_dir_prefix} ${custom_suffixes} ${remote_option} ${remove_suffix}
		"""
}







workflow remote_fastq_input {
	take:
		fastq_ch

	main:

		transfer_fastqs(fastq_ch.collect())
		res_ch = transfer_fastqs.out.fastqs.flatten()

	emit:
		fastqs = res_ch
}


workflow remote_bam_input {
	take:
		bam_ch
	main:
		transfer_bams(bam_ch)
		res_ch = transfer_bams.out.bamfiles.flatten()
	emit:
		bamfiles = res_ch
}


workflow fastq_input {
	take:
		fastq_ch
	
	main:
		prepare_fastqs(fastq_ch.collect(), (params.remote_input_dir != null || params.remote_input_dir))

		fastq_ch = prepare_fastqs.out.fastqs
			.flatten()
			.map { file -> 
				def sample = file.getParent().getName()
				return tuple(sample, file)
			}.groupTuple(sort: true)
			.map { classify_sample(it[0], it[1]) }

	emit:
		fastqs = fastq_ch
}


workflow bam_input {
	take:
		bam_ch
	main:

		if (params.remote_input_dir) {
			bam_ch = remote_bam_input(bam_ch.collect())
		}

		bam_ch = bam_ch
			.map { file ->
				def sample = file.name.replaceAll(bam_suffix_pattern, "").replaceAll(/\.$/, "")
				return tuple(sample, file)
			}
			.groupTuple(sort: true)
			.map { classify_sample(it[0], it[1]) }

		if (params.do_bam2fq_conversion) {
			bam2fq(bam_ch)
			bam_ch = bam2fq.out.reads
				.map { classify_sample(it[0].id, it[1]) }
		}
	emit:
		bamfiles = bam_ch
}
