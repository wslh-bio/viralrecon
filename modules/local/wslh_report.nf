process WSLH_REPORT {
    label 'process_medium'

    publishDir = [
        path: { "${params.outdir}/wslh_report" },
        mode: params.publish_dir_mode,
        pattern: 'wslh_report.csv'
    ]

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.20.0--py27h7835474_0' :
        'quay.io/hdc-workflows/pysam:latest' }"

    input:
    path csv_variants
    path bam_files 
    path pangolin_reports
    path nextclade_reports

    output:
    path '*.csv'       , emit: csv

    when: //may not need this section
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    """
    viralrecon_postrun.py \\
        -s $csv_variants \\
        -b $bam_files \\
        -p $pangolin_reports \\
        -n $nextclade_reports \\
        -wv ${workflow.manifest.version}
    """
}