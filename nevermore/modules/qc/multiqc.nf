process multiqc {
    
    input:
    path(reports)
	path(multiqc_config)
	val(stage)

    output:
    path("reports/${stage}.multiqc_report.html")

    script:
    def send_report = (false && params.email && params.mailer) ? "echo . | ${params.mailer} -s 'multiqc_report' -a reports/${stage}.multiqc_report.html ${params.email}" : ""
    """
	mkdir -p reports/
    multiqc -o reports/ -n ${stage}.multiqc_report.html -c ${multiqc_config} .
    ${send_report}
    """
}
