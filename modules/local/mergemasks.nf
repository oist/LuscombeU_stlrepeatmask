process MERGE_MASKS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(genome), path(tantan), path(windowmasker), path(repeatmasker)

    output:
    tuple val(meta), path("*.fasta.gz")    , emit: fasta
    tuple val(meta), path("*_jaccard.txt") , emit: txt
    tuple val(meta), path("*.bed.gz")      , emit: bed_gz
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools jaccard -nonamecheck -a $tantan       -b $windowmasker         > ${prefix}_tantan_windowmasker_jaccard.txt
    bedtools jaccard -nonamecheck -a $tantan       -b $repeatmasker         > ${prefix}_tantan_repeatmasker_jaccard.txt
    bedtools jaccard -nonamecheck -a $repeatmasker -b $windowmasker         > ${prefix}_repeatmasker_windowmasker_jaccard.txt

    zcat $tantan $windowmasker               | sort -k1,1 -k2,2n | bedtools merge | gzip --best --no-name > ${prefix}_tantan_windowmasker.bed.gz
    zcat $tantan $repeatmasker               | sort -k1,1 -k2,2n | bedtools merge | gzip --best --no-name > ${prefix}_tantan_repeatmasker.bed.gz
    zcat $windowmasker $repeatmasker         | sort -k1,1 -k2,2n | bedtools merge | gzip --best --no-name > ${prefix}_windowmasker_repeatmasker.bed.gz
    zcat $tantan $windowmasker $repeatmasker | sort -k1,1 -k2,2n | bedtools merge | gzip --best --no-name > ${prefix}_allmaskers.bed.gz

    bedtools \\
        maskfasta \\
        -soft \\
        -fi $genome \\
        -bed ${prefix}_allmaskers.bed.gz \\
        -fo /dev/stdout |
        gzip --best --no-name > ${prefix}_allmaskers.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
