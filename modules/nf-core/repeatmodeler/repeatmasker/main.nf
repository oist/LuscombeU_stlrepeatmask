process REPEATMODELER_MASKER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.5--pl5321hdfd78af_0':
        'biocontainers/repeatmodeler:2.0.5--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.fa") , emit: fastas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    RepeatMasker \\
        -lib $fasta \\
        $ref \\
        $args \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repeatmodeler: \$(RepeatModeler --version | sed 's/RepeatModeler version //')
    END_VERSIONS
    """
}
