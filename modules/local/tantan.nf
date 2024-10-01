process TANTAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tantan:49--h43eeafb_0':
        'biocontainers/tantan:49--h43eeafb_0' }"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path("*.masked.fa.gz"), emit: masked_fa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tantan ${ref} | gzip --best --no-name > ${prefix}.masked.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version |& sed '1!d ; s/tantan //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.masked.fa.gz
    gzip --no-name ${prefix}.masked.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version |& sed '1!d ; s/tantan //')
    END_VERSIONS
    """
}
