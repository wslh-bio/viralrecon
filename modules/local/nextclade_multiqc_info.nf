process NEXTCLADE_MULTIQC_INFO {
    tag "$db"
    label 'process_medium'

    conda "conda-forge::python=3.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8' :
        'quay.io/biocontainers/python:3.8' }"

    input:
    path db
    path clade_tsv

    output:
    path "nextclade_clade_db_info_mqc.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    nextclade_multiqc_info.py $db nextclade_clade_db_info_mqc.tsv $clade_tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
