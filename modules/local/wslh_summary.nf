process WSLH_SUMMARY {
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path csv_variants

    output:
    path("*_viralrecon_report.csv"), emit: wslh_report

    when:
    task.ext.when == null || task.ext.when

    script:
    def runname = params.runname ?: workflow.runName
    """
    wslh_report_summary.py \\
        -s ${csv_variants} \\
        -r ${runname} \\
        -wv ${workflow.manifest.version}
    """
}