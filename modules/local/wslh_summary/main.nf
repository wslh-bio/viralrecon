process WSLH_SUMMARY {
    label 'process_low'

    container "quay.io/wslh-bioinformatics/pandas@sha256:bf3cb8e5f695cc7c4cf8cc5ab7e7924d1fc4c40dfbe7cb907110e93a7bf6f101"

    input:
    path multiqc_csv_report
    val runname
    path pangolin_report
    path nextclade_report

    output:
    path("*_wslh_report.csv"), emit: qc_report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wslh_results_summary.py \\
    -qc ${multiqc_csv_report} \\
    -r ${runname} \\
    -p ${pangolin_report} \\
    -n ${nextclade_report} \\
    -wv ${workflow.version} \\
    """
}