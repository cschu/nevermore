process multiqc {
    input:
    path(reports)
	path(multiqc_config)
	val(stage)

    output:
    path("reports/${stage}.multiqc_report.html")

    script:
    """
	mkdir -p reports/
    multiqc -o reports/ -n ${stage}.multiqc_report.html -c ${multiqc_config} .
    """
}
