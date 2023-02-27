process remove_host_kraken2 {
	label 'kraken2'

    input:
    tuple val(sample), path(fq)
	path(kraken_db)

    output:
    tuple val(sample), path("no_host/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads

    script:
    def out_options = (sample.is_paired) ? "--paired --unclassified-out ${sample.id}#.fastq" : "--unclassified-out ${sample.id}_1.fastq"
    def move_r2 = (sample.is_paired) ? "gzip -c ${sample.id}_2.fastq > no_host/${sample.id}/${sample.id}_R2.fastq.gz" : ""

	def kraken2_call = "kraken2 --threads $task.cpus --db ${kraken_db} --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups}"

    """
    mkdir -p no_host/${sample.id}
	mkdir -p stats/decon/

    ${kraken2_call} ${out_options} --output stats/decon/${sample.id}.kraken_read_report.txt --report stats/decon/${sample.id}.kraken_report.txt $fq

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
	tuple val(sample), path("no_host/${sample.id}/${sample.id}_R*.fastq.gz"), emit: reads, optional: true
	tuple val(sample), path("no_host/${sample.id}/${sample.id}.chimeras_R1.fastq.gz"), emit: chimera_orphans, optional:true
	tuple val(sample), path("stats/decon/${sample.id}*.txt"), emit: stats, optional: true
	tuple val(sample), path("no_host/${sample.id}/KRAKEN_FINISHED"), emit: sentinel

	script:
	def kraken2_call = "kraken2 --threads $task.cpus --db ${kraken_db} --report-minimizer-data --gzip-compressed --minimum-hit-groups ${params.kraken2_min_hit_groups}"

	if (sample.is_paired) {
		"""
		set -e -o pipefail

		mkdir -p no_host/${sample.id}
		mkdir -p stats/decon/ 

		${kraken2_call} --unclassified-out ${sample.id}_1.fastq --output stats/decon/${sample.id}.kraken_read_report_1.txt --report stats/decon/${sample.id}.kraken_report_1.txt ${sample.id}_R1.fastq.gz
		${kraken2_call} --unclassified-out ${sample.id}_2.fastq --output stats/decon/${sample.id}.kraken_read_report_2.txt --report stats/decon/${sample.id}.kraken_report_2.txt ${sample.id}_R2.fastq.gz

		if [[ -f ${sample.id}_1.fastq || -f ${sample.id}_2.fastq ]]; then
			mkdir -p tmp/
			awk 'NR%4==1' *.fastq | sed 's/^@//' | cut -f 1 -d ' ' | sed 's/\\/[12]//' | sort -T tmp/ | uniq -c | sed 's/^\\s\\+//' > union.txt
			rm -rf tmp/

			((grep '^1' union.txt | cut -f 2 -d " ") || true) > chimeras.txt
			((grep '^2' union.txt | cut -f 2 -d " ") || true) > pairs.txt

			seqtk subseq ${sample.id}_1.fastq chimeras.txt | gzip -c - > chimeras.gz
			seqtk subseq ${sample.id}_1.fastq <(sed "s:\$:/1:" chimeras.txt) | gzip -c - >> chimeras.gz
			seqtk subseq ${sample.id}_2.fastq chimeras.txt | gzip -c - > chimeras.gz
			seqtk subseq ${sample.id}_2.fastq <(sed "s:\$:/2:" chimeras.txt) | gzip -c - >> chimeras.gz

			if [[ ! -z "\$(gzip -dc chimeras.gz | head -n 1)" ]]; then
				mv chimeras.gz no_host/${sample.id}/${sample.id}.chimeras_R1.fastq.gz
			fi

			seqtk subseq ${sample.id}_1.fastq pairs.txt | gzip -c - > no_host/${sample.id}/${sample.id}_R1.fastq.gz
			seqtk subseq ${sample.id}_1.fastq <(sed "s:\$:/1:" pairs.txt) | gzip -c - >> no_host/${sample.id}/${sample.id}_R1.fastq.gz
			seqtk subseq ${sample.id}_2.fastq pairs.txt | gzip -c - > no_host/${sample.id}/${sample.id}_R2.fastq.gz
			seqtk subseq ${sample.id}_2.fastq <(sed "s:\$:/2:" pairs.txt) | gzip -c - >> no_host/${sample.id}/${sample.id}_R2.fastq.gz

			rm -f ${sample.id}_1.fastq ${sample.id}_2.fastq
			rm -f chimeras.txt pairs.txt union.txt chimeras.gz
		fi

		touch no_host/${sample.id}/KRAKEN_FINISHED
		"""
	} else {
		"""
		mkdir -p no_host/${sample.id}
		mkdir -p stats/decon/

		${kraken2_call} --unclassified-out ${sample.id}_1.fastq --output stats/decon/${sample.id}.kraken_read_report_1.txt --report stats/decon/${sample.id}.kraken_report_1.txt ${sample.id}_R1.fastq.gz

		if [[ -f ${sample.id}_1.fastq ]]; then
			mv ${sample.id}_1.fastq no_host/${sample.id}/${sample.id}_R1.fastq
			gzip no_host/${sample.id}/*.fastq
		fi

		touch no_host/${sample.id}/KRAKEN_FINISHED
		"""
	}

}
