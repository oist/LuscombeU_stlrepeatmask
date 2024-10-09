process MERGE_REPM_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(genome), path(repeatmodeller), path(dfam), path(extlib)

    output:
    tuple val(meta), path("*.fasta.gz")                     , emit: fasta
    tuple val(meta), path("*_jaccard.txt")                  , emit: txt
    tuple val(meta), path("*.bed.gz")                       , emit: bed_gz
    tuple val(meta), path("*_repeatmasker_all.mask.bed.gz") , emit: repm_all_bed_gz
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '/^>/ {if (seqlen){print seqname "\t" seqlen}; split(\$1, a, ">"); seqname=a[2]; seqlen=0; next} {seqlen += length(\$0)} END {print seqname "\t" seqlen}' $genome > genome.genome # thanks, ChatGPT!
    run_bedtools_operations() {
        bedtools jaccard -nonamecheck -a "\$1" -b "\$2" -g genome.genome              > "${prefix}_\${3}_jaccard.txt"
        zcat "\$1" "\$2" | sort -k1,1 -k2,2n | bedtools merge | gzip --best --no-name > "${prefix}_\${3}.mask.bed.gz"
    }

    # Merge when available
    [ -e ${prefix}_DFAM.mask.bed.gz ]      && run_bedtools_operations "$repeatmodeller" "$dfam"   repeatmodeller_dfam
    [ -e ${prefix}_EXTR.mask.bed.gz ]      && run_bedtools_operations "$repeatmodeller" "$extlib" repeatmodeller_extr
    [ -e ${prefix}_EXTR.mask.bed.gz ] &&
        [ -e ${prefix}_DFAM.mask.bed.gz ]  && run_bedtools_operations "$dfam" "$extlib" dfam_extr
    [ -e ${prefix}_dfam_extr.mask.bed.gz ] && run_bedtools_operations "$repeatmodeller" ${prefix}_dfam_extr.mask.bed.gz all_combined

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
        gzip --best --no-name > ${prefix}_repeatmasker_all.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
