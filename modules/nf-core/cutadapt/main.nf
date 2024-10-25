process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.trim.fastq.gz" : "-o ${prefix}_1.trim.fastq.gz -p ${prefix}_2.trim.fastq.gz"
    def threeprime_adapters = params.threeprime_adapters ?: false

    // Ajustar primers según las condiciones dadas
    def primers
    if (params.threeprime_adapters && meta.single_end) {
        primers = "-a file:${meta.primers}"
    } else if (params.threeprime_adapters && !meta.single_end) {
        primers = "-a file:${meta.primers} -A file:${meta.primers}"
    } else if (!params.threeprime_adapters && meta.single_end) {
        primers = "-g file:${meta.primers}"
    } else {
        primers = "-g file:${meta.primers} -G file:${meta.primers}"
    }
    
    """
    cutadapt \\
        -Z \\
        --cores $task.cpus \\
        $args \\
        $primers \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trim.fastq.gz" : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    """
    touch ${prefix}.cutadapt.log
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
