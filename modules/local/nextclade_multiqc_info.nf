process NEXTCLADE_MULTIQC_INFO {
    tag "$db"
    label 'process_medium'

    conda "conda-forge::pandas:2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'quay.io/biocontainers/pandas:2.2.1' }"

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
    nextclade_multiqc_info.py --db_dir $db --clade_tsv $clade_tsv --out_tsv nextclade_clade_db_info_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
