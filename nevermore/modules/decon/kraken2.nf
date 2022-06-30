process remove_host_kraken2 {
	label 'kraken2'
	publishDir params.output_dir, mode: params.publish_mode

    input:
    tuple val(sample), path(fq)
	path(kraken_db)

    output:
    tuple val(sample), path("no_host/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

    script:
    def out_options = (sample.is_paired) ? "--paired --unclassified-out ${sample.id}#.fastq" : "--unclassified-out ${sample.id}_1.fastq"
    def move_r2 = (sample.is_paired) ? "gzip -c ${sample.id}_2.fastq > no_host/${sample.id}/${sample.id}_R2.fastq.gz" : ""

    """
    mkdir -p no_host/${sample.id}
    kraken2 --threads $task.cpus --db ${kraken_db} ${out_options} --output kraken_read_report.txt --report kraken_report.txt --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups} $fq

    gzip -c ${sample.id}_1.fastq > no_host/${sample.id}/${sample.id}_R1.fastq.gz
    ${move_r2}
    """
}


process remove_host_kraken2_individual {
	label 'kraken2'

	input:
	tuple val(sample), path(fq)
	path(kraken_db)

	output:
	tuple val(sample), path("no_host/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads
	tuple val(sample), path("no_host/${sample.id}/${sample.id}.c_orphans_R1.fastq.gz"), emit: chimera_orphans, optional:true

	script:
	def kraken2_call = "kraken2 --threads $task.cpus --db ${kraken_db} --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups}"

	if (sample.is_paired) {
		"""
		set -o pipefail

		mkdir -p no_host/${sample.id}
		mkdir -p ${sample.id}

		mkdir -p tmp/

		${kraken2_call} --unclassified-out ${sample.id}_1.fastq --output ${sample.id}/kraken_read_report_1.txt --report ${sample.id}/kraken_report_1.txt ${sample.id}_R1.fastq.gz
		${kraken2_call} --unclassified-out ${sample.id}_2.fastq --output ${sample.id}/kraken_read_report_2.txt --report ${sample.id}/kraken_report_2.txt ${sample.id}_R2.fastq.gz

		awk 'NR%4==1' *.fastq | sed 's/^@//' | cut -f 1 -d ' ' | sed 's/\\/[12]//' | sort -T tmp/ | uniq -c | sed 's/^\\s\\+//' > union.txt
		rm -rf tmp/
		grep '^1' union.txt | cut -f 2 -d " " > orphans.txt
		grep '^2' union.txt | cut -f 2 -d " " > pairs.txt

		seqtk subseq ${sample.id}_1.fastq orphans.txt | gzip -c - >> c_orphans.gz
		seqtk subseq ${sample.id}_2.fastq orphans.txt | gzip -c - >> c_orphans.gz

		if [[ ! -z "\$(gzip -dc c_orphans.gz | head -n 1)" ]]; then
			mv c_orphans.gz no_host/${sample.id}/${sample.id}.c_orphans_R1.fastq.gz
		fi

		seqtk subseq ${sample.id}_1.fastq pairs.txt | gzip -c - > no_host/${sample.id}/${sample.id}_R1.fastq.gz
		seqtk subseq ${sample.id}_2.fastq pairs.txt | gzip -c - > no_host/${sample.id}/${sample.id}_R2.fastq.gz
		"""
	} else {
		"""
		mkdir -p no_host/${sample.id}
		mkdir -p ${sample.id}

		${kraken2_call} --unclassified-out ${sample.id}_1.fastq --output ${sample.id}/kraken_read_report_1.txt --report ${sample.id}/kraken_report_1.txt ${sample.id}_R1.fastq.gz

		mv ${sample.id}_1.fastq no_host/${sample.id}/${sample.id}_R1.fastq
		gzip no_host/${sample.id}/*.fastq
		"""
	}

}
