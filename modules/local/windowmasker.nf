process WINDOWMASKER_MASK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.15.0--pl5321h6f7f691_1':
        'biocontainers/blast:2.15.0--pl5321h6f7f691_1' }"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.txt")          , emit: counts
    tuple val(meta), path("*.masked.fa.gz") , emit: masked_fa
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ""
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    windowmasker -mk_counts \\
        $args \\
        -mem ${(task.memory.toMega()).intValue()} \\
        -in ${ref} \\
        -out ${prefix}.txt

    windowmasker -ustat \\
        ${prefix}.txt \\
        $args \\
        -outfmt fasta \\
        -in ${ref} \\
        -out ${prefix}.masked.fa

    gzip --best --no-name ${prefix}.masked.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: "${meta.id}"

    """
    touch ${prefix}.masked.fa
    gzip --no-name ${prefix}.masked.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        windowmasker: \$(windowmasker -version-full | head -n 1 | sed 's/^.*windowmasker: //; s/ .*\$//')
    END_VERSIONS
    """
}
