process MERGE_REPM_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0d740f724375ad694bf4dce496aa7a419ffc67e12329bfb513935aafca5b28e9/data' :
        'community.wave.seqera.io/library/bedtools_blast_samtools_tantan:73b553483a4b3a4e' }"

    input:
    tuple val(meta), path(genome), path(repeatmodeller), path(extlib)

    output:
    tuple val(meta), path("*.fasta.gz")                     , emit: fasta
    tuple val(meta), path("*.bed.gz")                       , emit: bed_gz
    tuple val(meta), path("*_repeatmasker_all.mask.bed.gz") , emit: repm_all_bed_gz
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_bedtools_operations() {
        zcat "\$1" "\$2" | sort -k1,1 -k2,2n | bedtools merge | bgzip --threads $task.cpus --compress-level 9 > "${prefix}_\${3}.mask.bed.gz"
    }

    # Merge when available
    run_bedtools_operations "$repeatmodeller" "$extlib" repeatmodeller_extr

    # Pick the broadest combination available
    for file in all_combined repeatmodeller_extr repeatmodeller_dfam REPM; do
        [ -e "${prefix}_\${file}.mask.bed.gz" ] && cp "${prefix}_\${file}.mask.bed.gz" "${prefix}_repeatmasker_all.mask.bed.gz" && break
    done

    # Softmask the genome with the broadest combination
    bedtools \\
        maskfasta \\
        -soft \\
        -fi $genome \\
        -bed ${prefix}_repeatmasker_all.mask.bed.gz \\
        -fo /dev/stdout |
        bgzip --threads $task.cpus --compress-level 9 > ${prefix}_repeatmasker_all.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
