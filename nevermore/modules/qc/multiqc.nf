process multiqc {
    publishDir params.output_dir, mode: params.publish_mode

    input:
    path(reports)
	path(multiqc_config)

    output:
    path("multiqc_report.html")

    script:
    def send_report = (params.email && params.mailer) ? "echo . | ${params.mailer} -s 'multiqc_report' -a multiqc_report.html ${params.email}" : ""
    """
    multiqc -c ${multiqc_config} .
    ${send_report}
    """
}
