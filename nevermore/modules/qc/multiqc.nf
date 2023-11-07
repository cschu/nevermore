process multiqc {
    
    input:
    path(reports)
	path(multiqc_config)
	val(stage)

    output:
    path("reports/${stage}.multiqc_report.html")

    script:
    def send_report = (false && params.comms.email && params.comms.mailer) ? "echo . | ${params.comms.mailer} -s 'multiqc_report' -a reports/${stage}.multiqc_report.html ${params.comms.email}" : ""
    """
	mkdir -p reports/
    multiqc -o reports/ -n ${stage}.multiqc_report.html -c ${multiqc_config} .
    ${send_report}
    """
}
