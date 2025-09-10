process MAKE_VARIANTS_LONG_TABLE {

    conda "conda-forge::python=3.13.2 conda-forge::matplotlib=3.10.1 conda-forge::pandas=2.2.3 conda-forge::r-sys=3.4.3 conda-forge::regex=2024.11.6 conda-forge::scipy=1.15.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5ee6e81aff2205d76ad8755d2181f8ea1dd747daa92aaf9dcba943b69aa9f458/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_r-sys_pruned:23244d66110fcdf2' }"

    input:
    path bcftools_query, stageAs: "bcftools_query/*"
    path snpsift, stageAs: "snpsift/*"
    path pangolin, stageAs: "pangolin/*"

    output:
    path "*.csv"       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    """
    make_variants_long_table.py \\
        --bcftools_query_dir ./bcftools_query \\
        --snpsift_dir ./snpsift \\
        --pangolin_dir ./pangolin \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
